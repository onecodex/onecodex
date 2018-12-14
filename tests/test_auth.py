from requests.auth import HTTPBasicAuth

from onecodex import Api
from onecodex.lib.auth import BearerTokenAuth


def test_bearer_auth_from_env(api_data, monkeypatch):
    monkeypatch.setenv('ONE_CODEX_BEARER_TOKEN', 'mysecrettoken')
    ocx = Api(base_url='http://localhost:3000', cache_schema=True)
    assert isinstance(ocx._req_args['auth'], BearerTokenAuth)
    sample = ocx.Samples.get('761bc54b97f64980')
    assert sample.visibility != 'public'


def test_api_key_auth_from_env(api_data, monkeypatch):
    monkeypatch.setenv('ONE_CODEX_API_KEY', 'mysecretkey')
    ocx = Api(base_url='http://localhost:3000', cache_schema=True)
    assert isinstance(ocx._req_args['auth'], HTTPBasicAuth)


def test_bearer_auth_from_kwargs(api_data):
    ocx = Api(bearer_token='mysecrettoken', base_url='http://localhost:3000', cache_schema=True)
    assert isinstance(ocx._req_args['auth'], BearerTokenAuth)


def test_api_key_auth_from_kwargs(api_data):
    ocx = Api(api_key='mysecretkey', base_url='http://localhost:3000', cache_schema=True)
    assert isinstance(ocx._req_args['auth'], HTTPBasicAuth)
