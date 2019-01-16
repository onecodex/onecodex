from click.testing import CliRunner
import mock
import os
import pytest
from testfixtures import Replace

from onecodex import Cli
from tests.conftest import API_DATA

DATE_FORMAT = '%Y-%m-%d %H:%M'


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
def test_analysis_help(runner, api_data, mocked_creds_path):
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
    result = runner.invoke(Cli, ['classifications', '45a573fb7833449a', '--results'])
    assert result.exit_code == 0
    assert "Staphylococcus" in result.output


# Panels
def test_panel_instances(runner, api_data, mocked_creds_file):
    result = runner.invoke(Cli, ['panels'])
    assert result.exit_code == 0


# Samples
def test_samples(runner, api_data, mocked_creds_file):
    r0 = runner.invoke(Cli, ['samples'])
    r1 = runner.invoke(Cli, ['samples', '7428cca4a3a04a8e'])
    assert r0.exit_code == 0
    assert r1.exit_code == 0
    assert API_DATA['GET::api/v1/samples/7428cca4a3a04a8e']['$uri'] in r0.output
    assert API_DATA['GET::api/v1/samples/7428cca4a3a04a8e']['$uri'] in r1.output


# Login tests
def mock_fetch_api_key(username, password, server_url):
    return '123yuixha87yd87q3123uiqhsd8q2738'


def test_api_login(runner, mocked_creds_path):
    login_input = 'user@example.com' + '\n' + 'userpassword' + '\n'
    successful_login_msg = 'Your ~/.onecodex credentials file was successfully created.'
    with Replace('onecodex.auth.fetch_api_key_from_uname', mock_fetch_api_key):
        result = runner.invoke(Cli, ['login'], input=login_input)
        assert result.exit_code == 0
        assert successful_login_msg in result.output


@pytest.mark.parametrize('email,success,code', [
    ('incorrect+email@onecodex.com', False, 1),
    ('asmngs@onecodex.com', True, 0),
])
def test_api_key_login(runner, api_data, mocked_creds_path, email, success, code):
    result = runner.invoke(Cli, ['--api-key', '0' * 32, 'login'], input=email)
    assert result.exit_code == code
    assert os.path.exists(os.path.expanduser('~/.onecodex')) is success


def test_creds_file_exists(runner, mocked_creds_file):
    expected_message = "Credentials file already exists"

    result = runner.invoke(Cli, ["login"])
    assert result.exit_code == 0
    assert expected_message in result.output


def test_silent_login(runner, mocked_creds_file, api_data):
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
    result = runner.invoke(Cli, ['samples'], catch_exceptions=False)
    assert 'requires authentication' in result.output

    # bearer token in environment
    with mock.patch.dict(os.environ, {'ONE_CODEX_BEARER_TOKEN': '00000000000000000000000000000000'}):
        assert 'ONE_CODEX_BEARER_TOKEN' in os.environ
        assert 'ONE_CODEX_API_KEY' not in os.environ
        result = runner.invoke(Cli, ['samples'], catch_exceptions=False)
    assert 'requires authentication' not in result.output

    # bearer token in environment
    with mock.patch.dict(os.environ, {'ONE_CODEX_API_KEY': '00000000000000000000000000000000'}):
        assert 'ONE_CODEX_BEARER_TOKEN' not in os.environ
        assert 'ONE_CODEX_API_KEY' in os.environ
        result = runner.invoke(Cli, ['samples'], catch_exceptions=False)
    assert 'requires authentication' not in result.output


# Uploads
SEQUENCE = ('ACGTGTCGTAGGTAGCTACGACGTAGCTAACGTGTCGTAGCTACGACGTAGCTA'
            'ACGTGTCGTAGCTACGACGTAGCTAGGGACGTGTCGTAGCTACGACGTAGCTAG\n')


@pytest.mark.parametrize("files,threads", [
    (["temp.fa"], False),
    (["temp.fa"], True),
    (["temp.fa", "temp2.fa"], False),
    (["temp.fa", "temp2.fa"], True),
])
def test_standard_uploads(runner, mocked_creds_path, caplog, upload_mocks, files, threads):
    """Test single and multi file uploads, with and without threads
       (but not files >5GB)
    """
    with runner.isolated_filesystem():
        args = ['--api-key', '01234567890123456789012345678901', '-v', 'upload']
        if not threads:
            args += ['--max-threads', '1']
        for f in files:
            args.append(f)
            with open(f, mode='w') as f_out:
                f_out.write('>Test fasta\n')
                f_out.write(SEQUENCE)

        result = runner.invoke(Cli, args)
        assert result.exit_code == 0
        assert 'ab6276c673814123' in caplog.text  # mocked file id


def test_empty_upload(runner, mocked_creds_path, upload_mocks):
    with runner.isolated_filesystem():
        f = 'tmp.fa'
        f_out = open(f, mode='w')
        f_out.close()
        args = ['--api-key', '01234567890123456789012345678901', 'upload', f]
        result = runner.invoke(Cli, args)
        assert result.exit_code != 0


def test_paired_files(runner, mocked_creds_path, upload_mocks):
    with runner.isolated_filesystem():
        f, f2 = 'temp_R1.fa', 'temp_R2.fa'
        with open(f, mode='w') as f_out, open(f2, mode='w') as f_out2:
            f_out.write('>Test fasta\n')
            f_out.write(SEQUENCE)
            f_out2.write('>Test fasta\n')
            f_out2.write(SEQUENCE)

        args = ['--api-key', '01234567890123456789012345678901', 'upload', f, f2]
        # check that only one upload is kicked off for the pair of files
        patch1 = 'onecodex.lib.upload.upload_fileobj'
        patch2 = 'onecodex.lib.upload.FASTXInterleave'
        with mock.patch(patch1) as mp, mock.patch(patch2) as mp2:
            result = runner.invoke(Cli, args)
            assert mp.call_count == 1
            assert mp2.call_count == 1
        assert 'It appears there are paired files' in result.output
        assert result.exit_code == 0

        # Check with --forward and --reverse, should succeed
        args = ['--api-key', '01234567890123456789012345678901', 'upload', '--forward', f,
                '--reverse', f2]
        with mock.patch(patch1) as mp, mock.patch(patch2) as mp2:
            result4 = runner.invoke(Cli, args)
            assert mp.call_count == 1
            assert mp2.call_count == 1
        assert 'It appears there are paired files' not in result4.output
        assert result4.exit_code == 0

        # Check with only --forward, should fail
        args = ['--api-key', '01234567890123456789012345678901', 'upload', '--forward', f]
        result5 = runner.invoke(Cli, args)
        assert 'You must specify both forward and reverse files' in result5.output
        assert result5.exit_code != 0

        # Check with only --reverse, should fail
        args = ['--api-key', '01234567890123456789012345678901', 'upload', '--reverse', f2]
        result6 = runner.invoke(Cli, args)
        assert 'You must specify both forward and reverse files' in result6.output
        assert result6.exit_code != 0

        # Check with --forward, --reverse, and files, should fail
        args = ['--api-key', '01234567890123456789012345678901', 'upload', '--forward', f,
                '--reverse', f2, f2]
        with mock.patch(patch1) as mp, mock.patch(patch2) as mp2:
            result7 = runner.invoke(Cli, args)
        assert 'You may not pass a FILES argument' in result7.output
        assert result7.exit_code != 0
