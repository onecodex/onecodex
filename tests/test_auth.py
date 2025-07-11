import pytest
from requests.auth import HTTPBasicAuth

from onecodex import Api
from onecodex.lib.auth import BearerTokenAuth, check_version


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
