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


def make_creds_file():
    api_key = '123yuixha87yd87q3123uiqhsd8q2738'
    now = datetime.datetime.now().strftime(DATE_FORMAT)
    fake_creds = {'api_key': api_key, 'saved_at': now, 'updated_at': None}
    path = os.path.expanduser("~/.onecodex")
    with open(path, mode='w') as f:
        f.write(json.dumps(fake_creds))


def test_api_login(runner, mocked_creds_path):
    login_input = 'user@example.com' + '\n' + 'userpassword' + '\n'
    successful_login_msg = 'Your ~/.onecodex credentials file was successfully created.'
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


def test_logout_creds_dne(runner, mocked_creds_path):
    expected_message = "No One Codex API keys found."
    result = runner.invoke(Cli, ["logout"])
    assert result.exit_code == 1
    assert expected_message in result.output


# Uploads
SEQUENCE = ('ACGTGTCGTAGGTAGCTACGACGTAGCTAACGTGTCGTAGCTACGACGTAGCTA'
            'ACGTGTCGTAGCTACGACGTAGCTAGGGACGTGTCGTAGCTACGACGTAGCTAG\n')


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
            args += ['--max-threads', '1']
        for f in files:
            args.append(f)
            with open(f, mode='w') as f_out:
                f_out.write('>Test fasta\n')
                f_out.write(SEQUENCE)

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


def test_paired_files(runner, upload_mocks):
    import mock

    with runner.isolated_filesystem():
        f, f2 = 'temp_R1.fa', 'temp_R2.fa'
        with open(f, mode='w') as f_out, open(f2, mode='w') as f_out2:
            f_out.write('>Test fasta\n')
            f_out.write(SEQUENCE)
            f_out2.write('>Test fasta\n')
            f_out2.write(SEQUENCE)

        args = ['--api-key', '01234567890123456789012345678901', 'upload', f, f2]
        # check that only one upload is kicked off for the pair of files
        patch1 = 'onecodex.lib.upload.upload_file'
        patch2 = 'onecodex.lib.inline_validator.FASTXTranslator.close'
        with mock.patch(patch1) as mp, mock.patch(patch2) as mp2:
            result = runner.invoke(Cli, args)
            assert mp.call_count == 1
            assert mp2.call_count == 0  # We close in the upload_file call
        assert 'It appears there are paired files' in result.output
        assert result.exit_code == 0

        # Check with validate=False, should fail
        args = ['--api-key', '01234567890123456789012345678901', 'upload', f, f2,
                '--do-not-validate']
        result2 = runner.invoke(Cli, args)
        assert result2.exit_code != 0

        # Check with validate=False, interleave=False, should success
        args = ['--api-key', '01234567890123456789012345678901', 'upload', f, f2,
                '--do-not-validate', '--do-not-interleave']
        with mock.patch(patch1) as mp:
            result3 = runner.invoke(Cli, args)
            assert mp.call_count == 2
        assert result3.exit_code == 0


def test_large_uploads(runner, upload_mocks, monkeypatch):
    # a lot of funky mocking
    import mock

    def mockfilesize(path):
        if 'large' in path:
            return 5 * 1000 * 1000 * 1000 + 1
        else:
            return 500  # small

    monkeypatch.setattr(os.path, 'getsize', mockfilesize)

    with runner.isolated_filesystem():
        big_file = "large.fa"
        with open(big_file, mode='w') as f:
            f.write('>BIG!!!\n')
            f.write(SEQUENCE)

        args = ['--api-key', '01234567890123456789012345678901', 'upload', big_file]

        def side_effect(*args, **kwargs):
            """Side effect to ensure FASTXValidator gets properly read
            """
            args[0].read()
            return None

        with mock.patch('onecodex.lib.upload.upload_large_file') as mp:
            mp.side_effect = side_effect
            result = runner.invoke(Cli, args)
            assert mp.call_count == 1

        assert result.exit_code == 0
        assert 'All complete.' in result.output
