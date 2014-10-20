from nose.tools import raises
from onecodex.cli import OneCodexArgParser
import onecodex.api_v0


# Note: We parse argv[1:] by default, so no prog name is required
# throughout the input_args lists in these tests

@raises(SystemExit)
def test_cli_help_exit_0():
    parser = OneCodexArgParser()
    input_args = ["--help"]
    parser.parse_args(input_args)


@raises(SystemExit)
def test_cli_upload_req():
    # Check for required args
    parser = OneCodexArgParser()
    input_args = ["upload"]
    parser.parse_args(input_args)


def test_cli_upload():
    parser = OneCodexArgParser()
    input_args = ["upload", "my_file.fasta", "--no-threads"]
    args = parser.parse_args(input_args)
    assert isinstance(args.file, list)
    assert args.file[0] == "my_file.fasta"
    assert args.threads is False
    assert args.max_threads == onecodex.api_v0.DEFAULT_THREADS


def test_cli_samples():
    parser = OneCodexArgParser()
    input_args = ["samples", "uuid1", "uuid2", "uuid3"]
    args = parser.parse_args(input_args)
    assert isinstance(args.samples, list)
    assert len(args.samples) == 3


def test_cli_analyses():
    parser = OneCodexArgParser()
    input_args = ["analyses", "uuid1", "uuid2", "uuid3"]
    args = parser.parse_args(input_args)
    assert isinstance(args.analyses, list)
    assert len(args.analyses) == 3


@raises(SystemExit)
def test_cli_analyses_raw_error():
    parser = OneCodexArgParser()
    input_args = ["analyses", "uuid1", "uuid2", "uuid3", "--raw"]
    args = parser.parse_args(input_args)
    args.run(args)


@raises(SystemExit)
def test_cli_analyses_table_and_raw_error():
    parser = OneCodexArgParser()
    input_args = ["analyses", "uuid1", "--table", "--raw"]
    args = parser.parse_args(input_args)
    args.run(args)


@raises(SystemExit)
def test_cli_analyses_table_error():
    parser = OneCodexArgParser()
    input_args = ["analyses", "uuid1", "uuid2", "--table"]
    args = parser.parse_args(input_args)
    args.run(args)


def test_cli_references():
    parser = OneCodexArgParser()
    input_args = ["references"]
    args = parser.parse_args(input_args)
    assert isinstance(args.references, list)
    assert len(args.references) == 0
