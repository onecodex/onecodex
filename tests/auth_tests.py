import datetime
import mock
import os
from nose.tools import raises
from onecodex.auth import OneCodexAuth
from onecodex.cli import OneCodexArgParser


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
