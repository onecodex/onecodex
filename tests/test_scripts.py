# flake8: noqa
import hashlib
import os
import pytest

pytest.importorskip("skbio")
import shutil

from onecodex import Cli


@pytest.mark.parametrize(
    "paired,subset_pairs_independently,with_children,exclude_reads,include_lowconf",
    [
        # single end
        (False, False, False, False, False),
        (False, False, False, True, False),  # --exclude-reads
        (False, False, False, False, True),  # --include-lowconf
        (False, False, False, True, True),  # --exclude-reads --include-lowconf
        (False, False, True, False, False),  # --with-children
        (False, False, True, True, False),  # --with-children --exclude-reads
        (False, False, True, False, True),  # --with-children --include-lowconf
        (False, False, True, True, True),  # --with-children --exclude-reads --include-lowconf
        # paired-end
        (True, False, False, False, False),
        (True, False, False, True, False),  # --exclude-reads
        (True, False, False, False, True),  # --include-lowconf
        (True, False, False, True, True),  # --exclude-reads --include-lowconf
        (True, False, True, False, False),  # --with-children
        (True, False, True, True, False),  # --with-children --exclude-reads
        (True, False, True, False, True),  # --with-children --include-lowconf
        (True, False, True, True, True),  # --with-children --exclude-reads --include-lowconf
        (True, True, False, False, False),  # --subset-pairs-independently
        (True, True, False, True, False),  # --subset-pairs-independently --exclude-reads
        (True, True, False, False, True),  # --subset-pairs-independently --include-lowconf
        (
            True,
            True,
            False,
            True,
            True,
        ),  # --subset-pairs-independently --exclude-reads --include-lowconf
        (True, True, True, False, False),  # --subset-pairs-independently --with-children
        (
            True,
            True,
            True,
            True,
            False,
        ),  # --subset-pairs-independently --with-children --exclude-reads
        (
            True,
            True,
            True,
            False,
            True,
        ),  # --subset-pairs-independently --with-children --include-lowconf
        (
            True,
            True,
            True,
            True,
            True,
        ),  # --subset-pairs-independently --with-children --exclude-reads --include-lowconf
    ],
)
def test_subset_reads(
    runner,
    api_data,
    mocked_creds_file,
    paired,
    subset_pairs_independently,
    with_children,
    exclude_reads,
    include_lowconf,
):
    basedir = os.path.abspath(os.path.dirname(__file__))
    data_dir = os.path.join(basedir, "data/files")
    files = [
        "test_paired_filtering_001.fastq.gz.results.tsv.gz",
        "test_paired_filtering_R1_001.fastq.gz",
        "test_paired_filtering_R2_001.fastq.gz",
        "test_single_filtering_001.fastq.gz",
        "test_single_filtering_001.fastq.gz.results.tsv.gz",
    ]

    digests = {
        # single-end
        (False, False, False, False, False): ["cc4856f34c4e4742414de15b429466cd"],
        (False, False, False, True, False): ["df8b67bf7e4e92122387c0a46a2f0ad4"],
        (False, False, False, False, True): ["26fda71d0ce44292c87c0fb71222f3a9"],
        (False, False, False, True, True): ["bb2f55962087097a6ecea10e921fa137"],
        (False, False, True, False, False): ["8a773e43459561db491ec36a7d9e7df5"],
        (False, False, True, True, False): ["60473e16683b258cc7e4da7e733adb68"],
        (False, False, True, False, True): ["6081c9985bb9f1ccf32bdcf2b41da14d"],
        (False, False, True, True, True): ["bd009027ef68c8f0b3d5c31b38893ca7"],
        # paired-end
        (True, False, False, False, False): [
            "f483a5e69f32839f90aea13be55e2c73",
            "3e1fb9a6df1407d79733e40466fc143e",
        ],
        (True, False, False, True, False): [
            "478e2bbabf6ec09ddd9856da3176e238",
            "adc4e91104e4be9736a23b4b3a8d1b5e",
        ],
        (True, False, False, False, True): [
            "f483a5e69f32839f90aea13be55e2c73",
            "3e1fb9a6df1407d79733e40466fc143e",
        ],
        (True, False, False, True, True): [
            "a009992b18579e9f917478e22ba6d393",
            "ff8ea269b733ad0665ea75cd3bfd1b4c",
        ],
        (True, False, True, False, False): [
            "b4899aff26a94f4a5022252354cc32a3",
            "f7c97b5d094ad41db620374adbb9d9ba",
        ],
        (True, False, True, True, False): [
            "74811360d6203345e7f30eb6fbf34a92",
            "79c7956338fdb2431cdb28280ebf68e3",
        ],
        (True, False, True, False, True): [
            "b4899aff26a94f4a5022252354cc32a3",
            "f7c97b5d094ad41db620374adbb9d9ba",
        ],
        (True, False, True, True, True): [
            "3db7eca476d25fcc0e5adc37e8e83c56",
            "ceb7f837ab51d4bec8bb7a1d53f2ea86",
        ],
        (True, True, False, False, False): [
            "2a517f1f55d1e01a7dd2904121105af0",
            "ac1f2c8b35a8d269edc2a9ac3e1171c3",
        ],
        (True, True, False, True, False): [
            "252467b00875869305c30572c923e0a1",
            "f4081d90fc4ef90c23d8c53a847234dd",
        ],
        (True, True, False, False, True): [
            "2a517f1f55d1e01a7dd2904121105af0",
            "ac1f2c8b35a8d269edc2a9ac3e1171c3",
        ],
        (True, True, False, True, True): [
            "f548459a490b7579f8a989e2127caf2d",
            "cf0ea0fa04cf40b5120cfa5c3731e779",
        ],
        (True, True, True, False, False): [
            "225333af8f060124f43412661efc0dc1",
            "49df9e1764732d12aa4970720214b9bd",
        ],
        (True, True, True, True, False): [
            "74811360d6203345e7f30eb6fbf34a92",
            "79c7956338fdb2431cdb28280ebf68e3",
        ],
        (True, True, True, False, True): [
            "225333af8f060124f43412661efc0dc1",
            "49df9e1764732d12aa4970720214b9bd",
        ],
        (True, True, True, True, True): [
            "6411c6d6a1e3cceb43c55428d4e8cdb9",
            "3486c724c38b25a5bab32d7afac484bd",
        ],
    }

    def md5sum(filepath):
        checksum = hashlib.md5()
        with open(filepath, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                checksum.update(chunk)
        return checksum

    with runner.isolated_filesystem():
        for f in files:
            path = os.path.join(data_dir, f)
            dest = os.getcwd()
            shutil.copy(path, dest)

        args = ["scripts", "subset_reads"]

        if paired:
            args += [
                "bef0bc57dd7f4c43",
                "test_paired_filtering_R1_001.fastq.gz",
                "-r",
                "test_paired_filtering_R2_001.fastq.gz",
                "-t",
                "816",
            ]

            outfiles = [
                "test_paired_filtering_R1_001.filtered.fastq",
                "test_paired_filtering_R2_001.filtered.fastq",
            ]
        else:
            args += ["0f4ee4ecb3a3412f", "test_single_filtering_001.fastq.gz", "-t", "816"]

            outfiles = ["test_single_filtering_001.filtered.fastq"]

        if with_children:
            args += ["--with-children"]
            args[args.index("816")] = "2"

        if exclude_reads:
            args += ["--exclude-reads"]

        if subset_pairs_independently:
            args += ["--subset-pairs-independently"]

        if include_lowconf:
            args += ["--include-lowconf"]

        # run with validation
        result = runner.invoke(Cli, args, catch_exceptions=False)
        assert "Using cached read-level results" in result.output

        results_digests = []
        for f in outfiles:
            results_digests.append(md5sum(f).hexdigest())

        assert (
            results_digests
            == digests[
                (paired, subset_pairs_independently, with_children, exclude_reads, include_lowconf)
            ]
        )

        # run without validation (faster)
        result = runner.invoke(Cli, args + ["--do-not-validate"], catch_exceptions=False)
        assert "Using cached read-level results" in result.output

        results_digests = []
        for f in outfiles:
            results_digests.append(md5sum(f).hexdigest())

        assert (
            results_digests
            == digests[
                (paired, subset_pairs_independently, with_children, exclude_reads, include_lowconf)
            ]
        )
