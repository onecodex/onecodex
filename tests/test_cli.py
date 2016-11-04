# General imports
import datetime
import json
import os

# Testing imports
from click.testing import CliRunner
import pytest
from testfixtures import Replace

# Module imports
from tests.conftest import API_DATA
from onecodex import Cli

DATE_FORMAT = '%Y-%m-%d %H:%M'


@pytest.fixture(scope='function')
def runner():
    runner = CliRunner(env={"ONE_CODEX_API_BASE": "http://localhost:3000"})
    return runner


def test_cli_help(runner):
    for args in [None, '-h', '--help']:
        command = [args] if args is not None else None
        result = runner.invoke(Cli, command)
        assert result.exit_code == 0
        assert 'One Codex v1 API command line interface' in result.output


def test_version(runner):
    result = runner.invoke(Cli, ['--version'])
    assert result.exit_code == 0
    assert "onecodex, version" in result.output


# Test CLI without base override
def test_cli_wo_override(api_data, monkeypatch):
    monkeypatch.delattr("requests.sessions.Session.request")
    runner = CliRunner()
    result = runner.invoke(Cli, ['-v', '--help'])
    assert result.exit_code == 0
    assert 'One Codex v1 API command line interface' in result.output


# Analyses
def test_analysis_help(runner, api_data, mocked_creds_file):
    result = runner.invoke(Cli, ['analyses', '--help'])
    analysis_desc = "Retrieve performed analyses"
    assert result.exit_code == 0
    assert analysis_desc in result.output


def test_analyses(runner, api_data, mocked_creds_file):
    r0 = runner.invoke(Cli, ['analyses'])
    r1 = runner.invoke(Cli, ['analyses', '593601a797914cbf'])
    assert r0.exit_code == 0
    assert r1.exit_code == 0
    assert API_DATA['GET::api/v1/analyses/593601a797914cbf']['$uri'] in r0.output
    assert API_DATA['GET::api/v1/analyses/593601a797914cbf']['$uri'] in r1.output


# Classifications
def test_classification_instance(runner, api_data, mocked_creds_file):
    result = runner.invoke(Cli, ['classifications', '593601a797914cbf'])
    assert result.exit_code == 0
    assert API_DATA['GET::api/v1/classifications/593601a797914cbf']['$uri'] in result.output


def test_classifications_table(runner, api_data, mocked_creds_file, monkeypatch):
    result = runner.invoke(Cli, ['classifications', 'f9e4a5506b154953', '--results'])
    assert result.exit_code == 0
    assert "Salmonella" in result.output


# Panels
def test_panel_instances(runner, api_data, mocked_creds_file):
    result = runner.invoke(Cli, ['panels'])
    assert result.exit_code == 0


# Samples
def test_samples(runner, api_data, mocked_creds_file):
    r0 = runner.invoke(Cli, ['samples'])
    r1 = runner.invoke(Cli, ['samples', '761bc54b97f64980'])
    assert r0.exit_code == 0
    assert r1.exit_code == 0
    assert API_DATA['GET::api/v1/samples/761bc54b97f64980']['$uri'] in r0.output
    assert API_DATA['GET::api/v1/samples/761bc54b97f64980']['$uri'] in r0.output


# Login tests
def mock_fetch_api_key(username, password, server_url):
    return '123yuixha87yd87q3123uiqhsd8q2738'


def make_creds_file():
    api_key = '123yuixha87yd87q3123uiqhsd8q2738'
    now = datetime.datetime.now().strftime(DATE_FORMAT)
    fake_creds = {'api_key': api_key, 'saved_at': now, 'updated_at': None}
    path = os.path.expanduser("~/.onecodex")
    with open(path, mode='w') as f:
        f.write(json.dumps(fake_creds))


def test_api_login(runner, mocked_creds_file):
    login_input = 'user@example.com' + '\n' + 'userpassword' + '\n'
    successful_login_msg = 'Your ~/.onecodex credentials file successfully created.'
    with Replace('onecodex.auth.fetch_api_key_from_uname', mock_fetch_api_key):
        result = runner.invoke(Cli, ['login'], input=login_input)
        assert result.exit_code == 0
        assert successful_login_msg in result.output


def test_creds_file_exists(runner, mocked_creds_file):
    with runner.isolated_filesystem():
        make_creds_file()
        expected_message = "Credentials file already exists"

        result = runner.invoke(Cli, ["login"])
        assert result.exit_code == 0
        assert expected_message in result.output


def test_silent_login(runner, mocked_creds_file, api_data):
    with runner.isolated_filesystem():
        make_creds_file()
        result = runner.invoke(Cli, ['samples'])
        assert result.exit_code == 0


def test_creds_file_corrupted(runner, mocked_creds_file):
    path = os.path.expanduser("~/.onecodex")
    with open(path, mode='w') as f:
        f.write("aslkdjaslkd\nkasjdlkas\nasdkjaslkd908&S&&^")
    expected_message = "Your ~/.onecodex credentials file appears to be corrupted."

    result = runner.invoke(Cli, ["login"])
    assert result.exit_code == 1
    assert expected_message in result.output


def test_logout_creds_exists(runner, mocked_creds_file):
    with runner.isolated_filesystem():
        make_creds_file()
        expected_message = "Successfully removed One Codex credentials."
        path = os.path.expanduser("~/.onecodex")
        result = runner.invoke(Cli, ["logout"])
        assert result.exit_code == 0
        assert expected_message in result.output
        assert os.path.exists(path) is False


def test_logout_creds_dne(runner, mocked_creds_file):
    expected_message = "No One Codex API keys found."
    result = runner.invoke(Cli, ["logout"])
    assert result.exit_code == 1
    assert expected_message in result.output


# Uploads
@pytest.mark.parametrize("files,threads", [
    (["temp.fa"], False),
    (["temp.fa"], True),
    (["temp.fa", "temp2.fa"], False),
    (["temp.fa", "temp2.fa"], True),
])
def test_standard_uploads(runner, upload_mocks, files, threads):
    """Test single and multi file uploads, with and without threads
       (but not files >5GB)
    """
    with runner.isolated_filesystem():
        args = ['--api-key', '01234567890123456789012345678901', 'upload']
        if not threads:
            args.append('--no-threads')
        for f in files:
            args.append(f)
            with open(f, mode='w') as f_out:
                f_out.write('>Test fasta\n')
                f_out.write('ACGTAGCTAGCTGACTAGCTGACTGAC\n')

        result = runner.invoke(Cli, args)
        assert result.exit_code == 0
        assert 'ab6276c673814123' in result.output  # mocked file id


def test_empty_upload(runner, upload_mocks):
    with runner.isolated_filesystem():
        f = 'tmp.fa'
        f_out = open(f, mode='w')
        f_out.close()
        args = ['--api-key', '01234567890123456789012345678901', 'upload', f]
        result = runner.invoke(Cli, args)
        assert result.exit_code != 0


@pytest.mark.parametrize('awscli_installed,expected_code,expected_messages', [
    (False, 1, ['You must install the awscli package for files >5GB in size']),
    (True, 0, ['Starting large (>5GB) file upload. Please be patient while the file transfers...',
               'Successfully uploaded: large.fa']),
])
def test_large_uploads(runner, upload_mocks, monkeypatch,
                       awscli_installed, expected_code, expected_messages):
    # A bunch of imports for AWS CLI mocking
    import subprocess
    import sys
    from testfixtures.popen import MockPopen
    from mock import Mock

    def mockfilesize(path):
        if 'large' in path:
            return 5 * 1000 * 1000 * 1000 + 1
        else:
            return 500  # small

    monkeypatch.setattr(os.path, 'getsize', mockfilesize)

    if awscli_installed:
        mockpopen = MockPopen()
        cmd = ("AWS_ACCESS_KEY_ID=aws_key AWS_SECRET_ACCESS_KEY=aws_secret_key "
               "aws s3 cp large.fa s3://onecodex-multipart-uploads-encrypted/abcdef0987654321 "
               "--sse")
        mockpopen.set_command(cmd)
        monkeypatch.setattr(subprocess, 'Popen', mockpopen)
        sys.modules['awscli'] = Mock()

    with runner.isolated_filesystem():
        big_file = "large.fa"
        with open(big_file, mode='w') as f:
            f.write('>BIG!!!\n')
            f.write('ACGTGTCGTAGCTACGACGTAGCTAG\n')

        args = ['--api-key', '01234567890123456789012345678901', 'upload', big_file]
        result = runner.invoke(Cli, args)
        assert result.exit_code == expected_code
        for e in expected_messages:
            assert e in result.output
