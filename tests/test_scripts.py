import hashlib
import os
import shutil

import pytest

from onecodex import Cli
from tests.test_cli import make_creds_file


@pytest.mark.parametrize('paired,split_pairs', [
    (False, False),
    (True, True),
    (True, False),
])
def test_filter_reads(runner, api_data, mocked_creds_file, paired, split_pairs):
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
                '593601a797914cbf',
                'test_paired_filtering_R1_001.fastq.gz',
                '-r',
                'test_paired_filtering_R2_001.fastq.gz',
                '-t',
                '816'
            ]
            outfiles = ['test_paired_filtering_R1_001.filtered.fastq',
                        'test_paired_filtering_R2_001.filtered.fastq']
            digests = ['27f79512d9c3bb422b07e18e70963c82',
                       '14ed375841899c584ed42df686672b23']
            if split_pairs:
                args += ['--split-pairs']
                digests = ['5032a6287b0cae09d4cea032df92bddb',
                           '72fd9bff7774cc16c9ceccb66919fc28']
        else:
            args += [
                '5a4b7e3bd3a44006',
                'test_single_filtering_001.fastq.gz',
                '-t',
                '816',
            ]
            outfiles = ['test_single_filtering_001.filtered.fastq']
            digests = ['c8a2de041bc3025476d0cf2c566d926f']
        result = runner.invoke(Cli, args)

        assert 'Using cached read-level results' in result.output

        results_digests = []
        for f in outfiles:
            results_digests.append(md5sum(f).hexdigest())

        assert results_digests == digests
