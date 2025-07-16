from click.testing import CliRunner
import mock
import os
import os.path
import pytest
from testfixtures import Replace

from onecodex import Cli
from tests.conftest import API_DATA

DATE_FORMAT = "%Y-%m-%d %H:%M"


@pytest.fixture
def mock_file_upload():
    with mock.patch("onecodex.lib.upload._upload_sequence_fileobj") as m:
        yield m


@pytest.fixture
def mock_sample_get():
    with mock.patch("onecodex.models.sample.Samples.get") as m:
        yield m


# Click 8.2.0 changed exit code from 0 to 2 if no arguments are supplied
@pytest.mark.parametrize("args,exit_codes", [(None, {0, 2}), ("-h", {0}), ("--help", {0})])
def test_cli_help(runner, args, exit_codes):
    command = [args] if args is not None else None
    result = runner.invoke(Cli, command)
    assert result.exit_code in exit_codes
    assert "One Codex v1 API command line interface" in result.output


def test_version(runner):
    result = runner.invoke(Cli, ["--version"])
    assert result.exit_code == 0
    assert "onecodex, version" in result.output


# Test CLI without base override
def test_cli_wo_override(api_data, monkeypatch):
    monkeypatch.delattr("requests.sessions.Session.request")
    runner = CliRunner()
    result = runner.invoke(Cli, ["-v", "--help"])
    assert result.exit_code == 0
    assert "One Codex v1 API command line interface" in result.output


# Analyses
def test_analysis_help(runner, api_data, mocked_creds_path):
    result = runner.invoke(Cli, ["analyses", "--help"])
    analysis_desc = "Retrieve performed analyses"
    assert result.exit_code == 0
    assert analysis_desc in result.output


def test_analyses(runner, api_data, mocked_creds_file):
    r0 = runner.invoke(Cli, ["analyses"])
    r1 = runner.invoke(Cli, ["analyses", "593601a797914cbf"])
    assert r0.exit_code == 0
    assert r1.exit_code == 0
    assert API_DATA["GET::api/v1/analyses/593601a797914cbf"]["$uri"] in r0.output
    assert API_DATA["GET::api/v1/analyses/593601a797914cbf"]["$uri"] in r1.output


# Classifications
def test_classification_instance(runner, api_data, mocked_creds_file):
    result = runner.invoke(Cli, ["classifications", "593601a797914cbf"])
    assert result.exit_code == 0
    assert API_DATA["GET::api/v1/classifications/593601a797914cbf"]["$uri"] in result.output


def test_classifications_table(runner, api_data, mocked_creds_file, monkeypatch):
    result = runner.invoke(Cli, ["classifications", "45a573fb7833449a", "--results"])
    assert result.exit_code == 0
    assert "Staphylococcus" in result.output


# Documents
def test_documents_table(runner, api_data, mocked_creds_file):
    result = runner.invoke(Cli, ["documents", "list"])
    assert len(result.stdout.split("\n")) == 20
    assert "OneCodexTakeHome" in result.stdout

    result = runner.invoke(Cli, ["documents", "list", "--json"])
    assert len(result.stdout.split("\n")) == 258  # shorter because we ignore extra fields
    assert "OneCodexTakeHome" in result.stdout


# Panels
def test_panel_instances(runner, api_data, mocked_creds_file):
    result = runner.invoke(Cli, ["panels"])
    assert result.exit_code == 0


# Samples
def test_samples(runner, api_data, mocked_creds_file):
    r0 = runner.invoke(Cli, ["samples"])
    r1 = runner.invoke(Cli, ["samples", "7428cca4a3a04a8e"])
    assert r0.exit_code == 0
    assert r1.exit_code == 0
    assert API_DATA["GET::api/v1/samples/7428cca4a3a04a8e"]["$uri"] in r0.output
    assert API_DATA["GET::api/v1/samples/7428cca4a3a04a8e"]["$uri"] in r1.output


# Login tests
def mock_fetch_api_key(username, password, server_url):
    return "123yuixha87yd87q3123uiqhsd8q2738"


def test_api_login(runner, mocked_creds_path):
    login_input = "user@example.com" + "\n" + "userpassword" + "\n"
    successful_login_msg = "Your ~/.onecodex credentials file was successfully created."
    with Replace("onecodex.auth.fetch_api_key_from_uname", mock_fetch_api_key):
        result = runner.invoke(Cli, ["login"], input=login_input)
        assert result.exit_code == 0
        assert successful_login_msg in result.output


@pytest.mark.parametrize(
    "email,success,code",
    [("incorrect+email@onecodex.com", False, 1), ("asmngs@onecodex.com", True, 0)],
)
def test_api_key_login(runner, api_data, mocked_creds_path, email, success, code):
    result = runner.invoke(Cli, ["--api-key", "0" * 32, "login"], input=email)
    assert result.exit_code == code
    assert os.path.exists(os.path.expanduser("~/.onecodex")) is success


def test_creds_file_exists(runner, mocked_creds_file):
    expected_message = "Credentials file already exists"

    result = runner.invoke(Cli, ["login"])
    assert result.exit_code == 0
    assert expected_message in result.output


def test_silent_login(runner, mocked_creds_file, api_data):
    result = runner.invoke(Cli, ["samples"])
    assert result.exit_code == 0


def test_creds_file_corrupted(runner, mocked_creds_file):
    path = os.path.expanduser("~/.onecodex")
    with open(path, mode="w") as f:
        f.write("aslkdjaslkd\nkasjdlkas\nasdkjaslkd908&S&&^")
    expected_message = "Your ~/.onecodex credentials file appears to be corrupted."

    result = runner.invoke(Cli, ["login"])
    assert result.exit_code == 1
    assert expected_message in result.output


def test_logout_creds_exists(runner, mocked_creds_file):
    expected_message = "Successfully removed One Codex credentials."
    path = os.path.expanduser("~/.onecodex")
    result = runner.invoke(Cli, ["logout"])
    assert result.exit_code == 0
    assert expected_message in result.output
    assert os.path.exists(path) is False


def test_logout_creds_dne(runner, mocked_creds_path):
    expected_message = "No One Codex API keys found."
    result = runner.invoke(Cli, ["logout"])
    assert result.exit_code == 1
    assert expected_message in result.output


def test_auth_from_env(runner, api_data, mocked_creds_path):
    # no authentication method, no stored login
    result = runner.invoke(Cli, ["samples"], catch_exceptions=False)
    assert "requires authentication" in result.output

    # bearer token in environment
    with mock.patch.dict(
        os.environ, {"ONE_CODEX_BEARER_TOKEN": "00000000000000000000000000000000"}
    ):
        assert "ONE_CODEX_BEARER_TOKEN" in os.environ
        assert "ONE_CODEX_API_KEY" not in os.environ
        result = runner.invoke(Cli, ["samples"], catch_exceptions=False)
    assert "requires authentication" not in result.output

    # bearer token in environment
    with mock.patch.dict(os.environ, {"ONE_CODEX_API_KEY": "00000000000000000000000000000000"}):
        assert "ONE_CODEX_BEARER_TOKEN" not in os.environ
        assert "ONE_CODEX_API_KEY" in os.environ
        result = runner.invoke(Cli, ["samples"], catch_exceptions=False)
    assert "requires authentication" not in result.output


# Uploads
@pytest.mark.parametrize(
    "files,threads",
    [
        (["temp.fa"], False),
        (["temp.fa"], True),
        (["temp.fa", "temp2.fa"], False),
        (["temp.fa", "temp2.fa"], True),
    ],
)
def test_standard_uploads(
    runner, mocked_creds_path, caplog, generate_fastq, upload_mocks, files, threads
):
    """Test single and multi file uploads, with and without threads
    (but not files >5GB)
    """
    with (
        runner.isolated_filesystem(),
        mock.patch("onecodex.models.Projects.get", side_effect=lambda _: None),
        mock.patch("onecodex.lib.upload._s3_intermediate_upload") as s3_upload_mock,
    ):
        s3_upload_mock.return_value = {"sample_id": "7428cca4a3a04a8e"}
        args = [
            "--api-key",
            "01234567890123456789012345678901",
            "-v",
            "upload",
            "--project",
            "testproj",
            "--tag",
            "isolate",
            "--metadata",
            "totalige=1234",
        ]
        if not threads:
            args += ["--max-threads", "1"]
        for f in files:
            args.append(generate_fastq(f))

        result = runner.invoke(Cli, args, catch_exceptions=False, input="\n")
        assert result.exit_code == 0
        assert "7428cca4a3a04a8e" in caplog.text  # mocked file id


def test_empty_upload(runner, mocked_creds_path, upload_mocks):
    with runner.isolated_filesystem():
        f = "tmp.fa"
        f_out = open(f, mode="w")
        f_out.close()
        args = ["--api-key", "01234567890123456789012345678901", "upload", f]
        result = runner.invoke(Cli, args, catch_exceptions=False)
        assert result.exit_code != 0


@pytest.mark.parametrize(
    "files,n_samples_uploaded,n_files_uploaded,n_paired_files,n_multiline_groups",
    [
        # 1 files, 1 sample
        (["test.fq"], 1, 1, 0, 0),
        # 2 files, 1 sample
        (["test_R1.fq", "test_R2.fq"], 1, 2, 2, 0),
        (["test_1.fq", "test_2.fq"], 1, 2, 2, 0),
        # 2 files, 2 lines each, 1 sample
        (
            [
                "test_S1_L001_R1_001.fastq.gz",
                "test_S1_L001_R2_001.fastq.gz",
                "test_S1_L002_R1_001.fastq.gz",
                "test_S1_L002_R2_001.fastq.gz",
            ],
            1,
            2,
            4,
            1,
        ),
        (["dir_r1_test/test_R1.fq", "dir_r1_test/test_R2.fq"], 1, 2, 2, 0),
        (["test_1_1.fq", "test_1_2.fq"], 1, 2, 2, 0),
        # 3 files, 2 samples
        (["test_R1.fq", "test_R2.fq", "other.fq"], 2, 3, 2, 0),
    ],
)
def test_paired_and_multiline_files(
    runner,
    generate_fastq,
    mock_file_upload,
    mock_sample_get,
    mocked_creds_path,
    upload_mocks,
    files,
    n_samples_uploaded,
    n_files_uploaded,
    n_paired_files,
    n_multiline_groups,
):
    files = [generate_fastq(x) for x in files]

    args = ["--api-key", "01234567890123456789012345678901", "upload"] + files
    # check that 2 uploads are kicked off for the pair of files
    result = runner.invoke(Cli, args, catch_exceptions=False, input="\n\n\n")
    assert mock_file_upload.call_count == n_files_uploaded
    assert mock_sample_get.call_count == n_samples_uploaded
    assert result.exit_code == 0

    paired_files_prompt = "It appears there are {} paired files (of {} total)".format(
        n_paired_files, len(files)
    )
    if n_paired_files > 0:
        assert paired_files_prompt in result.output
    else:
        assert paired_files_prompt not in result.output

    multilane_prompt = "This data appears to have been split across multiple sequencing lanes.\nConcatenate lanes before upload?"
    if n_multiline_groups > 0:
        assert multilane_prompt in result.output
    else:
        assert multilane_prompt not in result.output


def test_paired_files_with_forward_and_reverse_args(
    runner, generate_fastq, mock_file_upload, mock_sample_get, mocked_creds_path, upload_mocks
):
    f1, f2 = generate_fastq("test_R1.fq"), generate_fastq("test_R2.fq")

    # Check with --forward and --reverse, should succeed
    args = [
        "--api-key",
        "01234567890123456789012345678901",
        "upload",
        "--forward",
        f1,
        "--reverse",
        f2,
    ]
    result = runner.invoke(Cli, args, input="Y")
    assert mock_file_upload.call_count == 2
    assert mock_sample_get.call_count == 1
    assert "It appears there are 2 paired files" not in result.output  # skips message
    assert result.exit_code == 0

    # Check with only --forward, should fail
    args = ["--api-key", "01234567890123456789012345678901", "upload", "--forward", f1]
    result = runner.invoke(Cli, args)
    assert "You must specify both forward and reverse files" in result.output
    assert result.exit_code != 0

    # Check with only --reverse, should fail
    args = ["--api-key", "01234567890123456789012345678901", "upload", "--reverse", f2]
    result = runner.invoke(Cli, args)
    assert "You must specify both forward and reverse files" in result.output
    assert result.exit_code != 0

    # Check with --forward, --reverse, and files, should fail
    args = [
        "--api-key",
        "01234567890123456789012345678901",
        "upload",
        "--forward",
        f1,
        "--reverse",
        f2,
        f2,
    ]
    result = runner.invoke(Cli, args)
    assert "You may not pass a FILES argument" in result.output
    assert result.exit_code != 0


@pytest.mark.parametrize(
    "files,n_samples_uploaded,n_files_uploaded,n_paired_files,n_ont_files",
    [
        # 1 files, 1 sample
        (["test.fq"], 1, 1, 0, 0),
        # 2 files, 1 sample
        (["test_1.fq", "test_2.fq"], 1, 2, 2, 0),
        # 4 files, 4 ONT parts, 1 sample
        (
            [
                "test_S1_0.fastq.gz",
                "test_S1_1.fastq.gz",
                "test_S1_2.fastq.gz",
                "test_S1_3.fastq.gz",
            ],
            1,
            1,
            0,
            4,
        ),
        # 3 files, 2 ONT parts, 2 samples
        (
            [
                "test_S1_0.fastq.gz",
                "test_S1_1.fastq.gz",
                "test_S2_0.fastq.gz",
            ],
            2,
            2,
            0,
            2,
        ),
        # 6 files, 2x3 ONT parts, 2 samples
        (
            [
                "test_S1_0.fastq.gz",
                "test_S2_0.fastq.gz",
                "test_S1_1.fastq.gz",
                "test_S2_1.fastq.gz",
                "test_S1_2.fastq.gz",
                "test_S2_2.fastq.gz",
            ],
            2,
            2,
            0,
            6,
        ),
        (["dir_r1_test/test_0.fq", "dir_r1_test/test_1.fq"], 1, 1, 0, 2),
        # 3 files, 2 samples
        (["test_0.fq", "other.fq", "test_1.fq"], 2, 2, 0, 2),
    ],
)
def test_paired_and_ont_files(
    runner,
    generate_fastq,
    mock_file_upload,
    mock_sample_get,
    mocked_creds_path,
    upload_mocks,
    files,
    n_samples_uploaded,
    n_files_uploaded,
    n_paired_files,
    n_ont_files,
):
    files = [generate_fastq(x) for x in files]

    args = ["--api-key", "01234567890123456789012345678901", "upload"] + files
    result = runner.invoke(Cli, args, catch_exceptions=False, input="\n\n")
    assert mock_file_upload.call_count == n_files_uploaded
    assert mock_sample_get.call_count == n_samples_uploaded
    assert result.exit_code == 0

    paired_files_prompt = "It appears there are {} paired files (of {} total)".format(
        n_paired_files, len(files)
    )
    if n_paired_files > 0:
        assert paired_files_prompt in result.output
    else:
        assert paired_files_prompt not in result.output

    ont_prompt = f"It appears there are {n_samples_uploaded} sample(s)"
    if n_ont_files > 0:
        assert ont_prompt in result.output
    else:
        assert ont_prompt not in result.output


def test_download_samples_without_prompt(runner, api_data, mocked_creds_file):
    with runner.isolated_filesystem():
        result = runner.invoke(Cli, ["download", "samples", "output", "--no-prompt"])
        assert result.exit_code == 0
        assert len(os.listdir("output")) > 0


def test_download_samples_with_prompt_confirmation(runner, api_data, mocked_creds_file):
    with runner.isolated_filesystem():
        result = runner.invoke(Cli, ["download", "samples", "output", "--prompt"], input="y\n")
        assert result.exit_code == 0
        assert len(os.listdir("output")) > 0


def test_download_samples_with_prompt_abort(runner, api_data, mocked_creds_file):
    with runner.isolated_filesystem():
        result = runner.invoke(Cli, ["download", "samples", "output", "--prompt"], input="n\n")
        assert result.exit_code == 0
        assert not os.path.exists("output")
