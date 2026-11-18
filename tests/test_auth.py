import datetime
import json

import mock
import pytest
import responses
from requests.auth import HTTPBasicAuth

from onecodex import Api
from onecodex.lib.auth import (
    DEVICE_NAME_HEADER,
    BearerTokenAuth,
    DeviceAuthorization,
    DeviceLoginError,
    check_version,
    poll_device_token,
    start_device_login,
    token_expires_at,
    token_time_left,
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


def test_expired_token_from_creds_file_is_ignored(api_data, mocked_token_creds_file, monkeypatch):
    mocked_token_creds_file(datetime.timedelta(minutes=-1))
    monkeypatch.setenv("ONE_CODEX_API_KEY", "mysecretkey")
    with pytest.warns(UserWarning, match="has expired"):
        ocx = Api(base_url="http://localhost:3000")
    assert isinstance(ocx._client.session.auth, HTTPBasicAuth)


@pytest.mark.parametrize(
    "creds,expected",
    [
        ({"api_key": "mysecretkey"}, None),
        (
            {"access_token": "t", "expires_at": "2026-11-23T15:00:00+00:00"},
            datetime.datetime(2026, 11, 23, 15, tzinfo=datetime.timezone.utc),
        ),
        # A timestamp without a timezone is UTC
        (
            {"access_token": "t", "expires_at": "2026-11-23T15:00:00"},
            datetime.datetime(2026, 11, 23, 15, tzinfo=datetime.timezone.utc),
        ),
    ],
)
def test_token_expires_at(creds, expected):
    assert token_expires_at(creds) == expected


@pytest.mark.parametrize("expires_at", [None, "", "next tuesday"])
def test_token_without_valid_expiry_is_expired(expires_at):
    creds = {"access_token": "t"}
    if expires_at is not None:
        creds["expires_at"] = expires_at
    assert token_time_left(creds) < datetime.timedelta(0)


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
