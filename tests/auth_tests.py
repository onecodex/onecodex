import datetime
import mock
import os
from nose.tools import raises
from onecodex.auth import OneCodexAuth
from onecodex.cli import OneCodexArgParser
from subprocess import check_output, CalledProcessError

TEST_CREDS = "onecodex_testing_creds.temp"


def test_auth_w_32_len_cli_api_key():
    parser = OneCodexArgParser()
    API_KEY = "12345678123456781234567812345678"
    input_args = ["upload", "my_file.fasta",
                  "--api-key", API_KEY]
    args = parser.parse_args(input_args)
    assert not hasattr(args, "credentials")
    OneCodexAuth(args, check_for_update=False, creds_file=TEST_CREDS)
    assert hasattr(args, "credentials")
    assert args.credentials["api_key"] == API_KEY
    assert args.credentials["saved_at"] is None
    assert args.credentials["updated_at"] is None


def test_auth_login():
    parser = OneCodexArgParser()
    API_KEY = "12345678123456781234567812345678"
    input_args = ["login"]
    with mock.patch("getpass.getpass", return_value=API_KEY):
        args = parser.parse_args(input_args)
        assert not hasattr(args, "credentials")
        OneCodexAuth(args, check_for_update=False, creds_file=TEST_CREDS)
    assert hasattr(args, "credentials")
    assert args.credentials["api_key"] == API_KEY
    assert args.credentials["saved_at"] == datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    assert args.credentials["updated_at"] is None  # We disable checking in the test


@raises(SystemExit)
def test_auth_login_file_already_exists():
    parser = OneCodexArgParser()
    input_args = ["login"]
    args = parser.parse_args(input_args)
    OneCodexAuth(args, check_for_update=False, creds_file=TEST_CREDS)


@raises(SystemExit)  # exit w/ 0
def test_auth_logout():
    parser = OneCodexArgParser()
    input_args = ["logout"]
    args = parser.parse_args(input_args)
    OneCodexAuth(args, check_for_update=False, creds_file=TEST_CREDS)
    assert not os.path.exists(TEST_CREDS)


@raises(SystemExit)
def test_auth_invalid_key_format():
    parser = OneCodexArgParser()
    input_args = ["upload", "my_file.fasta",
                  "--api-key",
                  "my_api_key_is_too_short"]
    args = parser.parse_args(input_args)
    OneCodexAuth(args, check_for_update=False, creds_file=TEST_CREDS)


# test aws cli proper failing with bad credentials
@raises(CalledProcessError)
def test_bad_aws_creds():
    
    # we modify the env here so we have to clean up
    replace_env = False
    aws_access = None
    aws_secret = None

    if os.environ.get("AWS_ACCESS_KEY_ID"):
        aws_access = os.environ["AWS_ACCESS_KEY_ID"]
        clean_env = True

    if os.environ.get("AWS_SECRET_ACCESS_KEY"):
        aws_secret = os.environ["AWS_SECRET_ACCESS_KEY"]
        clean_env = True

    # set bad credentials in ENV
    # then send the upload request
    try:
        os.environ["AWS_ACCESS_KEY_ID"] = "CRAPPYKEY" 
        os.environ["AWS_SECRET_ACCESS_KEY"] = "CRAPPYKEY"
        API_KEY = "12345678123456781234567812345678"
        parser = OneCodexArgParser()
        input_args = ["--api-key", API_KEY, "upload", "big1.fastq"]
        args = parser.parse_args(input_args)
        OneCodexAuth(args, check_for_update=False, creds_file=TEST_CREDS)
        assert hasattr(args, "credentials")
        assert args.credentials["api_key"] == API_KEY
        args.run(args)

    # catch the expected exit
    except SystemExit:

        # clean up the crappy env vars
        if replace_env:
            os.environ["AWS_ACCESS_KEY_ID"] = aws_access
            os.environ["AWS_SECRET_ACCESS_KEY"] = aws_secret
        else:
            del os.environ["AWS_ACCESS_KEY_ID"] 
            del os.environ["AWS_SECRET_ACCESS_KEY"]

        # try and find the pidof aws
        # aws should have been killed so this raises CalledProcessError
        check_output(["pidof","aws"])