from __future__ import print_function
from click.testing import CliRunner
from contextlib import contextmanager
import datetime
import gzip
import json
import mock
import os
from pkg_resources import resource_string
import pytest
import re
import requests
import responses

from onecodex import Api
from onecodex.lib.inline_validator import BaseFASTXReader

# we need to do this to make sure matplotlib isn't trying to connect to a GUI
try:
    import matplotlib
    matplotlib.use('Agg')
except ImportError:
    pass  # May not be incl. in minimal environment

# do not load EnhancedSampleCollection for tests unless explicitly needed
os.environ['ONE_CODEX_NO_ENHANCED'] = 'True'


def intercept(func, log=False, dump=None):
    """
    Used to copy API requests to make sure test data doesn't depend upon a connection to the One
    Codex server (basically like `betamax`, but for our requests/responses setup).

    For example, to dump out a log of everything that the function `test_function` requests, do the
    following:

    >>>mock_responses = {}
    >>>intercept(test_function, dump=mock_responses)
    >>>mock_json = json.dumps(mock_responses, separators=(',', ':'))

    Then you can test the function in the future by copying the output of mock_json into
    a string literal and doing:

    >>>mock_request(test_function, mock_json)
    """
    def handle_request(request):
        if log:
            print('->', request.method, request.url)

        # patch the request through (and disable mocking for this chunk)
        responses.mock.stop()
        resp = requests.get(request.url, headers=request.headers)
        text = resp.text
        headers = resp.headers
        # for some reason, responses pitches a fit about this being in the cookie
        headers['Set-Cookie'] = headers.get('Set-Cookie', '').replace(' HttpOnly;', '')
        responses.mock.start()
        data = json.dumps(json.loads(text), separators=(',', ':'))
        if log:
            print('<-', resp.status_code, data if log == 'all' else '*')
        if dump is not None:
            dump[request.method + ':' + request.url.split('/', 3)[-1]] = data
        return (200, headers, text)

    regex = re.compile('.*')
    with responses.mock as rsps:
        rsps.add_callback(responses.GET, regex, callback=handle_request)
        func()


# TODO: Fix a bug wherein this will return all the items to potion
#       but potion will still try to request subsequent pages... (it's stubborn!)
#       CRITICALLY THIS MEANS THAT TEST CASES SHOULD ONLY USE THE FIRST 20 ITEMS
@contextmanager
def mock_requests(mock_json):
    with responses.mock as rsps:
        for mock_url, mock_data in mock_json.items():
            method, content_type, url = mock_url.split(':', 2)
            if not content_type:
                content_type = 'application/json'
            if callable(mock_data):
                rsps.add_callback(method, re.compile('http://[^/]+/' + url + '(\\?.*)?$'),
                                  callback=mock_data,
                                  content_type=content_type)
            else:
                rsps.add(method, re.compile('http://[^/]+/' + url + '(\\?.*)?$'),
                         body=json.dumps(mock_data),
                         content_type=content_type)
        yield


def rs(path):
    return resource_string(__name__, path).decode('utf-8')


def json_resource(path):
    return json.loads(rs(path))


# All of the API data
# Scheme is
# METHOD:CONTENT_TYPE:URL  (content-type is optional)
# and then data is JSON or a callable
def _make_callback_resp(data, status_code=200):
    return status_code, {'Content-Type': 'application/json'}, json.dumps(data)


def make_potion_callback(status_code=200):
    def callback(req):
        data = json.loads(req.body)
        return _make_callback_resp(data, status_code=status_code)
    return callback


def update_metadata_callback(req):
    """
    Check that the metadata PATCH method isn't passing `sample`
    which will 400. Note this actually erases the sample attribute,
    which is kind of unfortunate, but a pain to patch without
    (we'd need to load the mocked metadata obj here and re-inject
    the sample)
    """
    metadata = json.loads(req.body)
    if 'sample' in metadata or '$uri' in metadata:  # write only field
        raise ValueError("Invalid metadata!")

    # Need to set $uri and sample for resource to resolve properly
    metadata['$uri'] = "/api/v1/metadata/4fe05e748b5a4f0e"
    metadata['sample'] = {"$ref": "/api/v1/samples/761bc54b97f64980"}
    return _make_callback_resp(metadata)


API_DATA = {
    # These are overrides for non-GET calls, which we don't auto-mock
    "DELETE::api/v1/samples/761bc54b97f64980": {},
    "GET::api/v1/samples/public": {},
    "GET::api/v1/projects/public": {},
    "GET::api/v1/analyses/public": [{
        "$uri": "/api/v1/analyses/de58fb3db76c42f3",
        "analysis_type": "classification",
        "complete": True,
        "created_at": "2017-08-03T17:23:21.931982+00:00",
        "error_msg": None,
        "job": {
            "$uri": "/api/v1/jobs/cc1d331e1ee54bac",
            "analysis_type": "classification",
            "created_at": "2016-05-05T17:27:02.116480+00:00",
            "name": "One Codex Database",
            "public": True
        },
        "sample": {
            "$uri": "/api/v1/samples/205bef56b3d6457d",
            "created_at": "2017-08-03T17:23:20.284145+00:00",
            "filename": "MSA-1000.16S.example.fastq.gz",
            "metadata": {
                "$ref": "/api/v1/metadata/ab296c81940e4cc0"
            },
            "owner": {
                "$ref": "/api/v1/users/5891ee65711c4d5e"
            },
            "primary_classification": {
                "$ref": "/api/v1/classifications/980c7388c54a46b8"
            },
            "project": {
                "$ref": "/api/v1/projects/66477c4c379a41a5"
            },
            "size": 9342972,
            "tags": [
                {
                    "$ref": "/api/v1/tags/e2910f93ec154617"
                },
                {
                    "$ref": "/api/v1/tags/7dd2d4ae9e9c4f14"
                },
                {
                    "$ref": "/api/v1/tags/7d350965f74346c1"
                },
                {
                    "$ref": "/api/v1/tags/1638c7a570214fda"
                }
            ],
            "visibility": "public"
        },
        "success": True,
    }],
    "GET::api/v1/classifications/45a573fb7833449a/results": {
        "table": [{
            'abundance': None,
            'name': 'Staphylococcus',
            'parent_tax_id': '1',
            'rank': 'genus',
            'readcount': 0,
            'readcount_w_children': 3,
            'tax_id': '1279'
        }, {
            'abundance': 1,
            'name': 'Staphylococcus sp. HGB0015',
            'parent_tax_id': '1279',
            'rank': 'species',
            'readcount': 3,
            'readcount_w_children': 3,
            'tax_id': '1078083',
        }]
    },
    "GET::api/v1/classifications/593601a797914cbf/results": {
        "table": [{
            'abundance': None,
            'name': 'Staphylococcus',
            'parent_tax_id': '1',
            'rank': 'genus',
            'readcount': 0,
            'readcount_w_children': 80,
            'tax_id': '1279'
        }, {
            'abundance': 1,
            'name': 'Staphylococcus sp. HGB0015',
            'parent_tax_id': '1279',
            'rank': 'species',
            'readcount': 80,
            'readcount_w_children': 80,
            'tax_id': '1078083',
        }]
    },
    "GET::api/v1/classifications/593601a797914cbf/readlevel": {
        "url": "https://s3.aws.com/bucket/test_paired_filtering_001.fastq.gz.results.tsv.gz"
    },
    "GET::api/v1/classifications/5a4b7e3bd3a44006/readlevel": {
        "url": "https://s3.aws.com/bucket/test_paired_filtering_001.fastq.gz.results.tsv.gz"
    },
    "GET::api/v1/classifications/bef0bc57dd7f4c43/readlevel": {
        "url": "https://s3.aws.com/bucket/test_paired_filtering_001.fastq.gz.results.tsv.gz"
    },
    "GET::api/v1/classifications/0f4ee4ecb3a3412f/readlevel": {
        "url": "https://s3.aws.com/bucket/test_single_filtering_001.fastq.gz.results.tsv.gz"
    },
    "PATCH::api/v1/samples/761bc54b97f64980": {},
    "PATCH::api/v1/metadata/4fe05e748b5a4f0e": update_metadata_callback
}

# explicitly load classification results for testing filter_reads
API_DATA['GET::api/v1/classifications/bef0bc57dd7f4c43/results'] = \
    json.load(open(os.path.join('tests', 'data', 'api', 'bef0bc57dd7f4c43_table.json')))
API_DATA['GET::api/v1/classifications/0f4ee4ecb3a3412f/results'] = \
    json.load(open(os.path.join('tests', 'data', 'api', '0f4ee4ecb3a3412f_table.json')))

# load classification results for testing viz
for c_id in ('6579e99943f84ad2', 'b50c176668234fe7', 'e0422602de41479f'):
    # converting bytes (returned by gzip) to string is necessary is python 3.4
    s = gzip.open(
        os.path.join('tests', 'data', 'api', '{}_table.json.gz'.format(c_id))
    ).read().decode('ascii')
    API_DATA['GET::api/v1/classifications/{}/results'.format(c_id)] = json.loads(s)

# and load the metadata that goes along with those samples
for md_id in ('3e7119ee74954abd', '6b69295478de4f1f', '6be1bb8849644f7b'):
    API_DATA['GET::api/v1/metadata/{}'.format(md_id)] = \
        json.load(open(os.path.join('tests', 'data', 'api', '{}_metadata.json'.format(md_id))))

for filename in os.listdir('tests/api_data'):
    if not filename.endswith('.json'):
        continue

    resource = json.load(open(os.path.join('tests/api_data', filename)))
    if filename.startswith('schema'):
        continue  # Parse separately below

    resource_name = filename.replace('.json', '')
    resource_uri = "GET::api/v1/{}".format(resource_name)
    API_DATA[resource_uri] = resource

    # Then iterate through all instances
    if isinstance(resource, list):
        for instance in resource:
            instance_uri = "GET::{}".format(instance['$uri'].lstrip('/'))
            API_DATA[instance_uri] = instance


SCHEMA_ROUTES = {}

for filename in os.listdir('tests/api_data'):
    if not filename.startswith('schema'):
        continue

    resource = json.load(open(os.path.join('tests/api_data', filename)))
    if filename == 'schema.json':
        resource_uri = 'GET::api/v1/schema'
    else:
        resource_name = filename.replace('.json', '').split('_')[1]
        resource_uri = 'GET::api/v1/{}/schema'.format(resource_name)

    SCHEMA_ROUTES[resource_uri] = resource


API_DATA.update(SCHEMA_ROUTES)


@pytest.yield_fixture(scope='function')
def api_data():
    with mock_requests(API_DATA):
        yield


@pytest.yield_fixture(scope='function')
def upload_mocks():
    def upload_callback(request):
        # Get and read the streaming iterator so it's empty
        if hasattr(request.body, 'fields'):
            streaming_iterator = request.body.fields['file'][1]
            streaming_iterator.read()
            assert isinstance(streaming_iterator, BaseFASTXReader)
        return (201, {'location': 'on-aws'}, '')

    json_data = {
        'POST:multipart/form-data:fake_aws_callback': upload_callback,
        'POST::api/v1/samples/init_upload': {
            'upload_url': 'http://localhost:3000/fake_aws_callback',
            'sample_id': 'ab6276c673814123',
            'additional_fields': {
                'AWSAccessKeyId': 'AKIAI36HUSHZTL3A7ORQ',
                'success_action_status': 201,
                'acl': 'private',
                'key': 'asd/file_ab6276c673814123/myfile.fastq',
                'signature': 'asdjsa',
                'policy': '123123123',
                'x-amz-server-side-encryption': 'AES256'
            }
        },
        'POST::api/v1/samples/confirm_upload': '',
        'GET::api/v1/samples/init_multipart_upload': {
            'callback_url': '/api/import_file_from_s3',
            'file_id': 'abcdef0987654321',
            's3_bucket': 'onecodex-multipart-uploads-encrypted',
            'upload_aws_access_key_id': 'aws_key',
            'upload_aws_secret_access_key': 'aws_secret_key'
        },
        'POST::api/import_file_from_s3': '',
    }
    json_data.update(SCHEMA_ROUTES)
    with mock_requests(json_data):
        yield


# API FIXTURES
@pytest.yield_fixture(scope='session')
def ocx_schemas():
    with mock_requests(SCHEMA_ROUTES):
        yield


@pytest.fixture(scope='session')
def ocx_w_raven():
    patched_env = os.environ.copy()
    patch = {
        'ONE_CODEX_API_BASE': 'http://localhost:3000',
        'ONE_CODEX_API_KEY': '1eab4217d30d42849dbde0cd1bb94e39',
        'ONE_CODEX_SENTRY_DSN': 'https://key:pass@sentry.example.com/1',
        'ONE_CODEX_NO_TELEMETRY': None,
    }

    patched_env.update(patch)

    with mock.patch.object(os, 'environ', patched_env):
        with mock_requests(SCHEMA_ROUTES):
            return Api(cache_schema=False, telemetry=True)


@pytest.fixture(scope='session')
def ocx():
    """Instantiated API client
    """
    with mock_requests(SCHEMA_ROUTES):
        return Api(api_key='1eab4217d30d42849dbde0cd1bb94e39',
                   base_url='http://localhost:3000', cache_schema=False)


@pytest.fixture(scope='function')
def runner():
    runner = CliRunner(env={
        'ONE_CODEX_API_BASE': 'http://localhost:3000',
        'ONE_CODEX_NO_TELEMETRY': 'True',
    })
    return runner


# CLI / FILE SYSTEM FIXTURE
@pytest.fixture(scope='function')
def mocked_creds_path(monkeypatch, tmpdir):
    # TODO: tmpdir is actually a LocalPath object
    # from py.path, and we coerce it into a string
    # for compatibility with the existing library code
    # *but* we should perhaps *not* do that for
    # better cross-platform compatibility. Investigate
    # and update as needed.
    def mockreturn(path):
        return os.path.join(str(tmpdir), '.onecodex')
    monkeypatch.setattr(os.path, 'expanduser', mockreturn)


@pytest.fixture(scope='function')
def mocked_creds_file(mocked_creds_path):
    with open(os.path.expanduser('~/.onecodex'), mode='w') as f:
        f.write(json.dumps({
            'email': 'test@onecodex.com',
            'api_key': None,
            'saved_at': datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
        }))
