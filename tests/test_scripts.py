import hashlib
import os
import pytest
import shutil

from onecodex import Cli
from tests.test_cli import make_creds_file


@pytest.mark.parametrize('paired,split_pairs,with_children,exclude_reads', [
    # single end
    (False, False, False, False),
    (False, False, True, False),  # --with-children
    (False, False, True, False),  # --exclude-reads
    (False, False, True, True),   # --with-children --exclude-reads

    # paired-end
    (True, False, False, False),
    (True, False, True, False),  # --with-children
    (True, False, False, True),  # --exclude-reads
    (True, False, True, True),   # --with-children --exclude-reads
    (True, True, False, False),  # --split-pairs
    (True, True, True, False),   # --split-pairs --with-children
    (True, True, False, True),   # --split-pairs --exclude-reads
    (True, True, True, True),    # --split-pairs --with-children --exclude-reads
])
def test_filter_reads(runner, api_data, mocked_creds_file, paired, split_pairs,
                      with_children, exclude_reads):
    make_creds_file()
    basedir = os.path.abspath(os.path.dirname(__file__))
    data_dir = os.path.join(basedir, 'data/files')
    files = [
        'test_paired_filtering_001.fastq.gz.results.tsv.gz',
        'test_paired_filtering_R1_001.fastq.gz',
        'test_paired_filtering_R2_001.fastq.gz',
        'test_single_filtering_001.fastq.gz',
        'test_single_filtering_001.fastq.gz.results.tsv.gz',
    ]

    def md5sum(filepath):
        checksum = hashlib.md5()
        with open(filepath, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b''):
                checksum.update(chunk)
        return checksum

    with runner.isolated_filesystem():
        for f in files:
            path = os.path.join(data_dir, f)
            dest = os.getcwd()
            shutil.copy(path, dest)

        args = ['scripts', 'filter_reads']
        if paired:
            args += [
                'bef0bc57dd7f4c43',
                'test_paired_filtering_R1_001.fastq.gz',
                '-r',
                'test_paired_filtering_R2_001.fastq.gz',
                '-t',
                '816'
            ]
            outfiles = ['test_paired_filtering_R1_001.filtered.fastq',
                        'test_paired_filtering_R2_001.filtered.fastq']
            digests = ['f483a5e69f32839f90aea13be55e2c73',
                       '3e1fb9a6df1407d79733e40466fc143e']
            if with_children:
                args += ['--with-children']
                args[args.index('816')] = '2'
                digests = ['b4899aff26a94f4a5022252354cc32a3',
                           'f7c97b5d094ad41db620374adbb9d9ba']

                if exclude_reads:
                    args += ['--exclude-reads']
                    digests = ['d6c512756d29f744e74c421c66ee95f1',
                               '33db92453df39be740ea1bb30c4e8be5']

                    if split_pairs:
                        args += ['--split-pairs']
                        digests = ['6411c6d6a1e3cceb43c55428d4e8cdb9',
                                   '3486c724c38b25a5bab32d7afac484bd']

            else:
                if exclude_reads:
                    args += ['--exclude-reads']
                    digests = ['f5aac827c45f24a138386ebbfc3e7611',
                               'a9fc38e273e766dcbc346445d6dac566']

                    if split_pairs:
                        args += ['--split-pairs']
                        digests = ['f548459a490b7579f8a989e2127caf2d',
                                   'cf0ea0fa04cf40b5120cfa5c3731e779']
        else:
            args += [
                '0f4ee4ecb3a3412f',
                'test_single_filtering_001.fastq.gz',
                '-t',
                '816',
            ]
            outfiles = ['test_single_filtering_001.filtered.fastq']
            digests = ['26fda71d0ce44292c87c0fb71222f3a9']
        result = runner.invoke(Cli, args)

        assert 'Using cached read-level results' in result.output

        results_digests = []
        for f in outfiles:
            results_digests.append(md5sum(f).hexdigest())

        assert results_digests == digests
