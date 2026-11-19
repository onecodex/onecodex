import datetime
import json
import os

import mock
import pytest
import requests
import responses
from requests.auth import HTTPBasicAuth

from onecodex import Api
from onecodex.lib.auth import (
    DEVICE_NAME_HEADER,
    BearerTokenAuth,
    DeviceAuthorization,
    DeviceLoginError,
    InvalidCredentialsError,
    OneCodexCredentials,
    check_version,
    poll_device_token,
    start_device_login,
)
from tests.conftest import DEVICE_AUTHORIZATION, DEVICE_TOKEN

SERVER = "http://localhost:3000/"


def test_bearer_auth_from_env(api_data, monkeypatch):
    monkeypatch.setenv("ONE_CODEX_BEARER_TOKEN", "mysecrettoken")
    ocx = Api(base_url="http://localhost:3000", cache_schema=True)
    assert isinstance(ocx._client.session.auth, BearerTokenAuth)
    sample = ocx.Samples.get("761bc54b97f64980")
    assert sample.visibility != "public"


def test_api_key_auth_from_env(api_data, monkeypatch):
    monkeypatch.setenv("ONE_CODEX_API_KEY", "mysecretkey")
    ocx = Api(base_url="http://localhost:3000", cache_schema=True)
    assert isinstance(ocx._client.session.auth, HTTPBasicAuth)


def test_bearer_auth_from_kwargs(api_data):
    ocx = Api(bearer_token="mysecrettoken", base_url="http://localhost:3000", cache_schema=True)
    assert isinstance(ocx._client.session.auth, BearerTokenAuth)


def test_api_key_auth_from_kwargs(api_data):
    ocx = Api(api_key="mysecretkey", base_url="http://localhost:3000", cache_schema=True)
    assert isinstance(ocx._client.session.auth, HTTPBasicAuth)


@pytest.mark.parametrize(
    "client_version,server_version,requires_update",
    [
        ("0.1.3", "0.1.3", False),
        ("0.1.4", "0.1.3", False),
        ("0.1.3", "0.1.4", True),
        ("0.2.3", "0.10.4", True),
    ],
)
def test_check_version(client_version, server_version, requires_update):
    from tests.conftest import mock_requests

    json_data = {"POST::api/v0/check_for_cli_update": {"latest_version": server_version}}

    server = "http://localhost:3000/"

    with mock_requests(json_data):
        result, error_code = check_version(client_version, server)
        assert result == requires_update


def test_token_auth_from_creds_file(api_data, mocked_token_creds_file):
    mocked_token_creds_file(datetime.timedelta(days=5))
    ocx = Api(base_url="http://localhost:3000")
    assert isinstance(ocx._client.session.auth, BearerTokenAuth)
    assert ocx._client.session.auth.token == "header.payload.signature"


def test_expiring_token_from_creds_file_warns(api_data, mocked_token_creds_file):
    mocked_token_creds_file(datetime.timedelta(hours=2))
    with pytest.warns(UserWarning, match="expires soon"):
        ocx = Api(base_url="http://localhost:3000")
    assert isinstance(ocx._client.session.auth, BearerTokenAuth)


def test_expired_token_from_creds_file_is_ignored(api_data, mocked_token_creds_file):
    mocked_token_creds_file(datetime.timedelta(minutes=-1))
    with pytest.warns(UserWarning, match="has expired"):
        ocx = Api(base_url="http://localhost:3000")
    assert ocx._client.session.auth is None


def test_api_key_auth_from_creds_file(api_data, mocked_creds_file):
    ocx = Api(base_url="http://localhost:3000")
    assert isinstance(ocx._client.session.auth, HTTPBasicAuth)
    assert ocx._client.session.auth.username == "123yuixha87yd87q3123uiqhsd8q2738"


@pytest.mark.parametrize("content", ["not json", "[]", '{"bearer_token": "t"}'])
def test_invalid_creds_file_warns(api_data, mocked_creds_path, content):
    with open(os.path.expanduser("~/.onecodex"), mode="w") as f:
        f.write(content)
    with pytest.warns(UserWarning, match="is corrupt"):
        ocx = Api(base_url="http://localhost:3000")
    assert ocx._client.session.auth is None


@pytest.mark.parametrize(
    "env,auth_type,secret",
    [
        ("ONE_CODEX_API_KEY", HTTPBasicAuth, "username"),
        ("ONE_CODEX_BEARER_TOKEN", BearerTokenAuth, "token"),
    ],
)
@pytest.mark.parametrize("token_creds", [True, False])
def test_env_takes_priority_over_creds_file(
    api_data,
    mocked_creds_file,
    mocked_token_creds_file,
    monkeypatch,
    env,
    auth_type,
    secret,
    token_creds,
):
    if token_creds:
        mocked_token_creds_file(datetime.timedelta(days=5))
    monkeypatch.setenv(env, "fromenv")
    ocx = Api(base_url="http://localhost:3000")
    assert isinstance(ocx._client.session.auth, auth_type)
    assert getattr(ocx._client.session.auth, secret) == "fromenv"


@pytest.mark.parametrize(
    "creds,expected",
    [
        ({"api_key": "mysecretkey"}, None),
        (
            {"bearer_token": "t", "expires_at": "2026-11-23T15:00:00+00:00"},
            datetime.datetime(2026, 11, 23, 15, tzinfo=datetime.timezone.utc),
        ),
        # A timestamp without a timezone is UTC
        (
            {"bearer_token": "t", "expires_at": "2026-11-23T15:00:00"},
            datetime.datetime(2026, 11, 23, 15, tzinfo=datetime.timezone.utc),
        ),
    ],
)
def test_credentials_expires_at(creds, expected):
    assert OneCodexCredentials(**creds).expires_at == expected


def test_credentials_time_left():
    expires_at = datetime.datetime.now(datetime.timezone.utc) + datetime.timedelta(hours=2)
    time_left = OneCodexCredentials(bearer_token="t", expires_at=expires_at).time_left
    assert datetime.timedelta(hours=1, minutes=59) < time_left <= datetime.timedelta(hours=2)
    assert OneCodexCredentials(api_key="mysecretkey").time_left is None


def test_credentials_save_and_load(tmp_path):
    path = str(tmp_path / ".onecodex")
    assert OneCodexCredentials.load(path) is None

    creds = OneCodexCredentials(
        email="test@onecodex.com",
        bearer_token="t",
        expires_at=datetime.datetime(2026, 11, 23, 15, tzinfo=datetime.timezone.utc),
        saved_at=datetime.datetime(2026, 11, 16, 15, 30),
    )
    creds.save(path)

    assert os.stat(path).st_mode & 0o777 == 0o600
    with open(path) as f:
        assert json.load(f) == {
            "email": "test@onecodex.com",
            "bearer_token": "t",
            "expires_at": "2026-11-23T15:00:00+00:00",
            "saved_at": "2026-11-16 15:30",
        }
    assert OneCodexCredentials.load(path) == creds


@pytest.mark.parametrize(
    "content",
    [
        "",
        "not json",
        "[]",
        # A login token requires a valid expiry
        '{"bearer_token": "t"}',
        '{"bearer_token": "t", "expires_at": ""}',
        '{"bearer_token": "t", "expires_at": "next tuesday"}',
        '{"api_key": 123}',
    ],
)
def test_credentials_load_invalid_file(tmp_path, content):
    path = tmp_path / ".onecodex"
    path.write_text(content)
    with pytest.raises(InvalidCredentialsError):
        OneCodexCredentials.load(str(path))


@pytest.mark.parametrize(
    "hostname,device_name",
    [
        ("MacBook-Pro-2.local", "MacBook-Pro-2"),
        ("lab-server-01", "lab-server-01"),
        ("Café", "Caf"),
        ("", None),
    ],
)
@responses.activate
def test_start_device_login(monkeypatch, hostname, device_name):
    monkeypatch.setattr("onecodex.lib.auth.socket.gethostname", lambda: hostname)
    responses.post(SERVER + "api/cli/device/code", json={**DEVICE_AUTHORIZATION, "extra": 1})

    assert start_device_login(SERVER) == DeviceAuthorization(**DEVICE_AUTHORIZATION)
    assert responses.calls[0].request.headers.get(DEVICE_NAME_HEADER) == device_name


@pytest.fixture
def mock_sleep():
    with mock.patch("onecodex.lib.auth.time.sleep") as mock_sleep:
        yield mock_sleep


@responses.activate
def test_poll_device_token(mock_sleep):
    url = SERVER + "api/cli/device/token"
    responses.post(url, status=400, json={"error": "authorization_pending"})
    responses.post(url, status=400, json={"error": "slow_down"})
    responses.post(url, json=DEVICE_TOKEN)

    token = poll_device_token(SERVER, DeviceAuthorization(**DEVICE_AUTHORIZATION))

    assert token.access_token == DEVICE_TOKEN["access_token"]
    assert token.email == DEVICE_TOKEN["email"]
    assert mock_sleep.call_args_list == [mock.call(5), mock.call(5), mock.call(10)]
    assert all(
        json.loads(call.request.body) == {"device_code": "secret-device-code"}
        for call in responses.calls
    )


@pytest.mark.parametrize(
    "timeout", [requests.exceptions.ReadTimeout, requests.exceptions.ConnectTimeout]
)
@responses.activate
def test_poll_device_token_retries_after_timeout(mock_sleep, timeout):
    url = SERVER + "api/cli/device/token"
    responses.post(url, body=timeout("No response in 15 seconds"))
    responses.post(url, json=DEVICE_TOKEN)

    token = poll_device_token(SERVER, DeviceAuthorization(**DEVICE_AUTHORIZATION))

    assert token.access_token == DEVICE_TOKEN["access_token"]
    # The wait doubles after the timeout
    assert mock_sleep.call_args_list == [mock.call(5), mock.call(10)]


@pytest.mark.parametrize(
    "error,message", [("access_denied", "denied"), ("expired_token", "expired")]
)
@responses.activate
def test_poll_device_token_failures(mock_sleep, error, message):
    responses.post(SERVER + "api/cli/device/token", status=400, json={"error": error})
    with pytest.raises(DeviceLoginError, match=message):
        poll_device_token(SERVER, DeviceAuthorization(**DEVICE_AUTHORIZATION))


@responses.activate
def test_poll_device_token_gives_up_after_expiry(mock_sleep, monkeypatch):
    responses.post(
        SERVER + "api/cli/device/token", status=400, json={"error": "authorization_pending"}
    )
    # Deadline at 600s, polls at 250s and 500s, then gives up at 750s
    clock = iter(range(0, 10_000, 250))
    monkeypatch.setattr("onecodex.lib.auth.time.monotonic", lambda: next(clock))

    with pytest.raises(DeviceLoginError, match="expired"):
        poll_device_token(SERVER, DeviceAuthorization(**DEVICE_AUTHORIZATION))
    assert len(responses.calls) == 2
