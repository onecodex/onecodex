import gzip
import os

import click
import pytest

from onecodex.input_helpers import (
    PlannedSample,
    confirm_plan,
    materialize_plan,
    plan_uploads,
)
from onecodex.utils import use_tempdir
from tests.conftest import FASTQ_SEQUENCE


def _names(samples):
    """Render a plan as (forward, reverse) basenames, for comparison in tests."""
    return [
        (
            tuple(os.path.basename(f) for f in sample.forward),
            tuple(os.path.basename(f) for f in sample.reverse) if sample.is_paired else None,
        )
        for sample in samples
    ]


@pytest.mark.parametrize(
    "files,expected",
    [
        # nothing to assemble
        (["test.fq"], [(("test.fq",), None)]),
        (["a.fq", "b.fq"], [(("a.fq",), None), (("b.fq",), None)]),
        # ONT chunks
        (["test_0.fq", "test_1.fq"], [(("test_0.fq", "test_1.fq"), None)]),
        (
            ["test_0.fq", "test_1.fq", "test_2.fq"],
            [(("test_0.fq", "test_1.fq", "test_2.fq"), None)],
        ),
        (["dir/test_0.fq", "dir/test_1.fq"], [(("test_0.fq", "test_1.fq"), None)]),
        # a chunk sequence has to start at 0 and have no gaps
        (["test_1.fq", "test_2.fq"], [(("test_1.fq",), ("test_2.fq",))]),
        (
            ["test_0.fq", "test_1.fq", "test_3.fq"],
            [(("test_0.fq",), None), (("test_1.fq",), None), (("test_3.fq",), None)],
        ),
        # a single chunk is a whole file already
        (["test_0.fq"], [(("test_0.fq",), None)]),
        # paired end reads
        (["test_R1.fq", "test_R2.fq"], [(("test_R1.fq",), ("test_R2.fq",))]),
        (["test_r1.fq", "test_r2.fq"], [(("test_r1.fq",), ("test_r2.fq",))]),
        (["test.R1.fq.gz", "test.R2.fq.gz"], [(("test.R1.fq.gz",), ("test.R2.fq.gz",))]),
        # one sample of each kind at once
        (
            ["s_0.fq", "s_1.fq", "p_R1.fq", "p_R2.fq", "other.fq"],
            [
                (("s_0.fq", "s_1.fq"), None),
                (("p_R1.fq",), ("p_R2.fq",)),
                (("other.fq",), None),
            ],
        ),
        # lanes
        (
            ["m_L001.fq", "m_L002.fq"],
            [(("m_L001.fq", "m_L002.fq"), None)],
        ),
        (
            ["m_L001_R1.fq", "m_L001_R2.fq", "m_L002_R1.fq", "m_L002_R2.fq"],
            [(("m_L001_R1.fq", "m_L002_R1.fq"), ("m_L001_R2.fq", "m_L002_R2.fq"))],
        ),
        # lanes have to start at 1 and have no gaps
        (
            ["m_L001.fq", "m_L003.fq"],
            [(("m_L001.fq",), None), (("m_L003.fq",), None)],
        ),
    ],
)
def test_plan_uploads(generate_fastq, files, expected):
    files = [generate_fastq(x) for x in files]
    assert sorted(_names(plan_uploads(files, prompt=False))) == sorted(expected)


def test_plan_uploads_finds_the_rest_of_a_chunk_sequence(generate_fastq):
    """Chunks that were not passed in are picked up from the same directory."""
    for filename in ["test_0.fq", "test_2.fq", "test_3.fq"]:
        generate_fastq(filename)
    files = [generate_fastq("test_1.fq")]

    assert _names(plan_uploads(files, prompt=True)) == [
        (("test_0.fq", "test_1.fq", "test_2.fq", "test_3.fq"), None)
    ]


def test_plan_uploads_finds_a_paired_mate(generate_fastq):
    generate_fastq("test_R2.fq")
    files = [generate_fastq("test_R1.fq")]

    assert _names(plan_uploads(files, prompt=True)) == [(("test_R1.fq",), ("test_R2.fq",))]


def test_plan_uploads_does_not_find_files_without_a_prompt(generate_fastq):
    """Files that were not passed in are only picked up if we can ask about them first."""
    for filename in ["test_0.fq", "test_R2.fq"]:
        generate_fastq(filename)
    files = [generate_fastq(x) for x in ["test_1.fq", "test_R1.fq"]]

    assert sorted(_names(plan_uploads(files, prompt=False))) == sorted(
        [(("test_1.fq",), None), (("test_R1.fq",), None)]
    )


def test_plan_uploads_ignores_a_file_cut_off_from_the_sequence(generate_fastq):
    """A chunk the run cannot reach must not drag in the chunks before the gap."""
    for filename in ["test_0.fq", "test_1.fq"]:
        generate_fastq(filename)
    files = [generate_fastq("test_5.fq")]

    assert _names(plan_uploads(files, prompt=True)) == [(("test_5.fq",), None)]


def test_plan_uploads_ignores_directories(generate_fastq, tmp_path):
    """A directory named like a chunk must not be treated as one."""
    files = [generate_fastq(x) for x in ["test_1.fq", "test_2.fq"]]
    os.mkdir(os.path.join(os.path.dirname(files[0]), "test_0.fq"))

    assert _names(plan_uploads(files, prompt=True)) == [(("test_1.fq",), ("test_2.fq",))]


def test_plan_uploads_separates_directories(generate_fastq):
    """Samples that share a filename but live in different directories stay separate."""
    files = [
        generate_fastq(x) for x in ["a/test_0.fq", "a/test_1.fq", "b/test_0.fq", "b/test_1.fq"]
    ]
    samples = plan_uploads(files, prompt=False)

    assert len(samples) == 2
    assert {tuple(os.path.dirname(f) for f in s.forward) for s in samples} == {
        (os.path.dirname(files[0]),) * 2,
        (os.path.dirname(files[2]),) * 2,
    }


def test_plan_uploads_never_uses_a_file_twice(generate_fastq):
    filenames = ["s_0.fq", "s_1.fq", "s_2.fq", "s_3.fq", "p_R1.fq", "p_R2.fq", "other.fq"]
    files = [generate_fastq(x) for x in filenames]

    used = [f for sample in plan_uploads(files, prompt=True) for f in sample.files]
    assert sorted(used) == sorted(files)


def test_materialize_plan_concatenates_in_order(generate_fastq):
    files = [generate_fastq(x) for x in ["test_0.fq", "test_1.fq", "test_2.fq"]]
    samples = plan_uploads(files, prompt=False)

    with use_tempdir() as tempdir:
        uploads = materialize_plan(samples, tempdir)
        assert len(uploads) == 1
        assert os.path.basename(uploads[0]) == "test.fq"
        with open(uploads[0]) as fin:
            assert fin.read() == 3 * FASTQ_SEQUENCE


def test_materialize_plan_concatenates_lanes_of_a_pair(generate_fastq_gz):
    files = [
        generate_fastq_gz(x)
        for x in ["m_L001_R1.fq.gz", "m_L001_R2.fq.gz", "m_L002_R1.fq.gz", "m_L002_R2.fq.gz"]
    ]
    samples = plan_uploads(files, prompt=False)

    with use_tempdir() as tempdir:
        ((forward, reverse),) = materialize_plan(samples, tempdir)
        assert os.path.basename(forward) == "m_R1.fq.gz"
        assert os.path.basename(reverse) == "m_R2.fq.gz"
        for path in (forward, reverse):
            with gzip.open(path, "r") as fin:
                assert fin.read() == 2 * FASTQ_SEQUENCE.encode("utf-8")


def test_materialize_plan_keeps_samples_with_the_same_name_apart(generate_fastq):
    files = [
        generate_fastq(x) for x in ["a/test_0.fq", "a/test_1.fq", "b/test_0.fq", "b/test_1.fq"]
    ]
    samples = plan_uploads(files, prompt=False)

    with use_tempdir() as tempdir:
        uploads = materialize_plan(samples, tempdir)
        assert [os.path.basename(u) for u in uploads] == ["test.fq", "test.fq"]
        assert len(set(uploads)) == 2


def test_materialize_plan_leaves_single_files_alone(generate_fastq):
    files = [generate_fastq("test.fq")]
    with use_tempdir() as tempdir:
        assert materialize_plan(plan_uploads(files, prompt=False), tempdir) == files


def test_confirm_plan_declined(generate_fastq, monkeypatch):
    monkeypatch.setattr(click, "prompt", lambda *args, **kwargs: "n")
    files = [generate_fastq(x) for x in ["test_R1.fq", "test_R2.fq"]]
    assert confirm_plan(plan_uploads(files, prompt=True), set(files)) is False


def test_confirm_plan_canceled(generate_fastq, monkeypatch):
    monkeypatch.setattr(click, "prompt", lambda *args, **kwargs: "c")
    files = [generate_fastq(x) for x in ["test_R1.fq", "test_R2.fq"]]
    with pytest.raises(SystemExit) as excinfo:
        confirm_plan(plan_uploads(files, prompt=True), set(files))
    assert excinfo.value.code == 0


def test_describe_plan_marks_files_found_on_disk(generate_fastq, monkeypatch, capsys):
    monkeypatch.setattr(click, "prompt", lambda *args, **kwargs: "Y")
    mate = generate_fastq("test_R2.fq")
    files = [generate_fastq("test_R1.fq")]

    assert confirm_plan(plan_uploads(files, prompt=True), set(files)) is True
    out = capsys.readouterr().out
    assert f"{mate} *" in out
    assert "* not specified on the command line" in out
    assert "1 sample(s) from 2 file(s)" in out


def test_planned_sample_shape():
    single = PlannedSample(forward=("a.fq",))
    assert not single.is_paired and not single.is_concatenated
    assert single.files == ("a.fq",)

    pair = PlannedSample(forward=("a_R1.fq",), reverse=("a_R2.fq",))
    assert pair.is_paired and not pair.is_concatenated
    assert pair.files == ("a_R1.fq", "a_R2.fq")
