import os
import gzip
import pytest
from onecodex.input_helpers import (
    _find_multiline_groups,
    concatenate_multiline_files,
    auto_detect_pairs,
)
from tests.conftest import FASTQ_SEQUENCE


def _get_basenames(elems):
    return [
        (os.path.basename(elem[0]), os.path.basename(elem[1]))
        if isinstance(elem, tuple)
        else os.path.basename(elem)
        for elem in elems
    ]


@pytest.mark.parametrize(
    "files,expected_pairing",
    [
        (["test.fq"], ["test.fq"]),
        (["test_R1.fq", "test_R2.fq"], [("test_R1.fq", "test_R2.fq")]),
        (["dir_r1_test/test_R1.fq", "dir_r1_test/test_R2.fq"], [("test_R1.fq", "test_R2.fq")]),
        (["test_1.fq", "test_2.fq"], [("test_1.fq", "test_2.fq")]),
        (["test_1_1.fq", "test_1_2.fq"], [("test_1_1.fq", "test_1_2.fq")]),
        (
            ["test_S1_L001_R1_001.fastq.gz", "test_S1_L001_R2_001.fastq.gz"],
            [("test_S1_L001_R1_001.fastq.gz", "test_S1_L001_R2_001.fastq.gz")],
        ),
        (["test_R1.fq", "test_R2.fq", "other.fq"], [("test_R1.fq", "test_R2.fq"), "other.fq"]),
    ],
)
def test_auto_detect_pairs(generate_fastq, files, expected_pairing):
    files = [generate_fastq(x) for x in files]
    pairs = auto_detect_pairs(files, prompt=False)
    basenames = _get_basenames(pairs)
    assert basenames == expected_pairing


def test_find_multiline_groups():
    files = [
        ("Sample1_L001_R1.fq", "Sample1_L001_R2.fq"),
        ("Sample1_L002_R1.fq", "Sample1_L002_R2.fq"),
        ("Sample1_L003_R1.fq", "Sample1_L003_R2.fq"),
        ("Sample2_L001_R1.fq", "Sample2_L001_R2.fq"),
        ("Sample2_L003_R1.fq", "Sample2_L003_R2.fq"),  # proper paired group
        ("Sample3_R1.fq", "Sample3_R2.fq"),  # no multiline
        ("Sample4_L001_R1.fq", "Sample4_L001_R2.fq"),
        "Sample4_L002_R2.fq",  # mismatch: R2 has more files
        "Sample5.fq",
        "Sample6_L001.fq",
        "Sample6_L002.fq",
        "Sample6_L003.fq",  # proper single group
        "Sample7_L001.fq",
        "Sample7_L002.fq",
        "Sample7_L004.fq",  # sequence gap
    ]
    expected_groups = [
        ["Sample6_L001.fq", "Sample6_L002.fq", "Sample6_L003.fq"],
        [
            ("Sample1_L001_R1.fq", "Sample1_L001_R2.fq"),
            ("Sample1_L002_R1.fq", "Sample1_L002_R2.fq"),
            ("Sample1_L003_R1.fq", "Sample1_L003_R2.fq"),
        ],
    ]

    groups = _find_multiline_groups(files)
    assert groups == expected_groups


def test_concatenate_multiline_files(generate_fastq):
    pairs = [
        (generate_fastq("Sample1_L001_R1.fq"), generate_fastq("Sample1_L001_R2.fq")),
        (generate_fastq("Sample1_L002_R1.fq"), generate_fastq("Sample1_L002_R2.fq")),
    ]
    singles = [
        generate_fastq("Sample2_L001.fq"),
        generate_fastq("Sample2_L002.fq"),
        generate_fastq("Sample2_L003.fq"),
    ]
    non_multiline = [("Sample3_R1.fq", "Sample3_R2.fq"), "Sample3.fq"]
    files = [*pairs, *singles, *non_multiline]

    concatenated = concatenate_multiline_files(files, prompt=False)

    basenames = _get_basenames(concatenated)
    assert basenames == [*non_multiline, "Sample2.fq", ("Sample1_R1.fq", "Sample1_R2.fq")]

    with open(concatenated[len(non_multiline)], "r") as inf:
        assert inf.read() == len(singles) * FASTQ_SEQUENCE

    with open(concatenated[len(non_multiline) + 1][0], "r") as inf:
        assert inf.read() == len(pairs) * FASTQ_SEQUENCE


def test_concatenate_gzipped_multiline_files(generate_fastq_gz):
    files = [
        generate_fastq_gz("Sample2_L001.fq"),
        generate_fastq_gz("Sample2_L002.fq"),
        generate_fastq_gz("Sample2_L003.fq"),
    ]
    concatenated = concatenate_multiline_files(files, prompt=False)
    assert len(concatenated) == 1
    with gzip.open(concatenated[0], "r") as fin:
        assert fin.read() == len(files) * FASTQ_SEQUENCE.encode("utf-8")
