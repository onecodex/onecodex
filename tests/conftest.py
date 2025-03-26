from __future__ import print_function
from click.testing import CliRunner
from contextlib import contextmanager
import copy
import datetime
import gzip
import json
import os
import pytest
import re
import responses

from onecodex import Api
from onecodex.analyses import FunctionalAnnotations, FunctionalAnnotationsMetric
from onecodex.models.collection import SampleCollection


# TODO: Fix a bug wherein this will return all the items to potion
#       but potion will still try to request subsequent pages... (it's stubborn!)
#       CRITICALLY THIS MEANS THAT TEST CASES SHOULD ONLY USE THE FIRST 20 ITEMS
@contextmanager
def mock_requests(mock_json):
    with responses.mock as rsps:
        for mock_url, mock_data in mock_json.items():
            method, content_type, url = mock_url.split(":", 2)
            if not content_type:
                content_type = "application/json"
            if "?" in url:
                compiled_url = re.compile("http://[^/]+/" + url)
            else:
                compiled_url = re.compile("http://[^/]+/" + url + "(\\?.*)?$")
            if callable(mock_data):
                rsps.add_callback(
                    method, compiled_url, callback=mock_data, content_type=content_type
                )
            else:
                rsps.add(
                    method, compiled_url, body=json.dumps(mock_data), content_type=content_type
                )
        yield rsps


def _make_callback_resp(data, status_code=200):
    return status_code, {"Content-Type": "application/json"}, json.dumps(data)


def update_metadata_callback(req):
    """
    Check that the metadata PATCH method isn't passing `sample`
    which will 400. Note this actually erases the sample attribute,
    which is kind of unfortunate, but a pain to patch without
    (we'd need to load the mocked metadata obj here and re-inject
    the sample)
    """
    metadata = json.loads(req.body)
    if "sample" in metadata or "$uri" in metadata:  # write only field
        raise ValueError("Invalid metadata!")

    # Need to set $uri and sample for resource to resolve properly
    metadata["$uri"] = "/api/v1/metadata/4fe05e748b5a4f0e"
    metadata["sample"] = {"$ref": "/api/v1/samples/761bc54b97f64980"}
    return _make_callback_resp(metadata)


def filtered_raw_results(raw_results, annotation, metric, taxa_stratified):
    metric = "total_" + metric if not taxa_stratified and annotation != "pathways" else metric
    table = [
        {
            "id": elem["id"] + "_" + elem["taxon_name"] if taxa_stratified else elem["id"],
            "name": elem["name"],
            "value": elem["value"],
            "taxon_id": elem["taxon_id"],
            "taxon_name": elem["taxon_name"],
        }
        for elem in raw_results["table"]
        if elem["group_name"] == annotation
        and elem["metric"] == metric
        and elem["taxa_stratified"] == taxa_stratified
    ]
    return {"table": table, "n_reads": raw_results["n_reads"], "n_mapped": raw_results["n_mapped"]}


# All of the mocked API data. Scheme is METHOD:CONTENT_TYPE:URL (content-type is optional) and then
# data is JSON or a callable. Things are also auto-mocked from the contents of tests/data/api/. See
# README.md for a description of how this all works.
API_DATA = {
    # These are overrides for non-GET calls, which we don't auto-mock
    "DELETE::api/v1/samples/761bc54b97f64980": {},
    "GET::api/v1/samples/public": [],
    "GET::api/v1/samples/organization": [
        {
            "$uri": "/api/v1/samples/7428cca4a3a04a8e",
            "created_at": "2015-09-25T17:27:19.596555-07:00",
            "filename": "SRR2352185.fastq.gz",
            "metadata": {"$ref": "/api/v1/metadata/a7fc7e430e704e2e"},
            "owner": {"$ref": "/api/v1/users/4ada56103d9a48b8"},
            "primary_classification": {"$ref": "/api/v1/classifications/464a7ebcf9f84050"},
            "project": None,
            "size": 181687821,
            "tags": [
                {"$ref": "/api/v1/tags/42997b7a62634985"},
                {"$ref": "/api/v1/tags/fb8e3b693c874f9e"},
                {"$ref": "/api/v1/tags/ff4e81909a4348d9"},
            ],
            "visibility": "private",
        }
    ],
    "GET::api/v1/projects/public": {},
    "GET::api/v1/analyses/public": [
        {
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
                "public": True,
            },
            "sample": {
                "$uri": "/api/v1/samples/205bef56b3d6457d",
                "created_at": "2017-08-03T17:23:20.284145+00:00",
                "filename": "MSA-1000.16S.example.fastq.gz",
                "metadata": {"$ref": "/api/v1/metadata/ab296c81940e4cc0"},
                "owner": {"$ref": "/api/v1/users/5891ee65711c4d5e"},
                "primary_classification": {"$ref": "/api/v1/classifications/980c7388c54a46b8"},
                "project": {"$ref": "/api/v1/projects/66477c4c379a41a5"},
                "size": 9342972,
                "tags": [
                    {"$ref": "/api/v1/tags/e2910f93ec154617"},
                    {"$ref": "/api/v1/tags/7dd2d4ae9e9c4f14"},
                    {"$ref": "/api/v1/tags/7d350965f74346c1"},
                    {"$ref": "/api/v1/tags/1638c7a570214fda"},
                ],
                "visibility": "public",
            },
            "success": True,
        }
    ],
    "GET::api/v1/classifications/45a573fb7833449a/results": {
        "table": [
            {
                "abundance": None,
                "name": "root",
                "parent_tax_id": None,
                "rank": "no rank",
                "readcount": 0,
                "readcount_w_children": 3,
                "tax_id": "1",
            },
            {
                "abundance": None,
                "name": "Staphylococcus",
                "parent_tax_id": "1",
                "rank": "genus",
                "readcount": 0,
                "readcount_w_children": 3,
                "tax_id": "1279",
            },
            {
                "abundance": 1,
                "name": "Staphylococcus sp. HGB0015",
                "parent_tax_id": "1279",
                "rank": "species",
                "readcount": 3,
                "readcount_w_children": 3,
                "tax_id": "1078083",
            },
        ]
    },
    "GET::api/v1/classifications/593601a797914cbf/results": {
        "table": [
            {
                "abundance": None,
                "name": "root",
                "parent_tax_id": None,
                "rank": "no rank",
                "readcount": 0,
                "readcount_w_children": 3,
                "tax_id": "1",
            },
            {
                "abundance": None,
                "name": "Staphylococcus",
                "parent_tax_id": "1",
                "rank": "genus",
                "readcount": 0,
                "readcount_w_children": 80,
                "tax_id": "1279",
            },
            {
                "abundance": 1,
                "name": "Staphylococcus sp. HGB0015",
                "parent_tax_id": "1279",
                "rank": "species",
                "readcount": 80,
                "readcount_w_children": 80,
                "tax_id": "1078083",
            },
        ]
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
    "PATCH::api/v1/metadata/4fe05e748b5a4f0e": update_metadata_callback,
    "POST::api/v1/samples/.*/download_uri": {
        "download_uri": "http://localhost:3000/mock/download/url"
    },
    "GET::mock/download/url": "1234567890",
    "POST::api/v1/documents/a4f6727a840a4df0/download_uri": {
        "download_uri": "http://localhost:3000/mock/download/url"
    },
    "GET::api/v1/jobs/cc1d331e1ee54bac": {
        "$uri": "/api/v1/jobs/cc1d331e1ee54bac",
        "analysis_type": "classification",
        "created_at": "2016-05-05T17:27:02.116480+00:00",
        "name": "One Codex Database (2017)",
        "public": True,
    },
    "GET::api/v1_experimental/jobs\\?.*where=%7B%22%24uri%22%3A\\+%7B%22%24in%22%3A\\+%5B%22%2Fapi%2Fv1_experimental%2Fjobs%2F59e7904ea8ed4202%22%5D%7D%7D&sort=%7B%22created_at%22%3A\\+true%7D": [
        {
            "$uri": "/api/v1_experimental/jobs/59e7904ea8ed4202",
            "analysis_type": "functional",
            "created_at": "2023-04-28T15:27:40.140791-07:00",
            "name": "Functional v1",
            "public": True,
        }
    ],
    "GET::api/v1/projects/4b53797444f846c4": {
        "$uri": "/api/v1/projects/472fc57510e24150",
        "description": None,
        "name": "Test",
        "owner": {"$ref": "/api/v1/users/9923090af03c46ce"},
        "permissions": [
            "can_see_files",
            "can_incur_charges",
            "can_download_files",
            "can_edit_metadata",
            "can_add_files",
            "can_administer",
        ],
        "project_name": "testproj",
        "public": False,
    },
    "GET::api/v1_experimental/projects/4b53797444f846c4": {
        "$uri": "/api/v1_experimental/projects/472fc57510e24150",
        "description": None,
        "name": "Test",
        "owner": {"$ref": "/api/v1_experimental/users/9923090af03c46ce"},
        "permissions": [
            "can_see_files",
            "can_incur_charges",
            "can_download_files",
            "can_edit_metadata",
            "can_add_files",
            "can_administer",
        ],
        "project_name": "testproj",
        "public": False,
    },
    "GET::api/v1/samples\\?.*where=%7B%22project%22%3A\\+%224b53797444f846c4%22%7D.*": [
        {
            "$uri": "/api/v1/samples/0b2d0b5397324841",
            "created_at": "2017-06-23T23:52:51.201676+00:00",
            "filename": "SRR4408293.fastq",
            "metadata": {"$ref": "/api/v1/metadata/6b69295478de4f1f"},
            "owner": {"$ref": "/api/v1/users/7189f36afe3640ac"},
            "primary_classification": {"$ref": "/api/v1/classifications/6579e99943f84ad2"},
            "project": {"$ref": "/api/v1/projects/4b53797444f846c4"},
            "size": 2225829803,
            "tags": [],
            "visibility": "public",
        },
        {
            "$uri": "/api/v1/samples/1640864a28bf44ba",
            "created_at": "2017-06-23T23:56:19.556902+00:00",
            "filename": "SRR4305031.fastq",
            "metadata": {"$ref": "/api/v1/metadata/6be1bb8849644f7b"},
            "owner": {"$ref": "/api/v1/users/7189f36afe3640ac"},
            "primary_classification": {"$ref": "/api/v1/classifications/b50c176668234fe7"},
            "project": {"$ref": "/api/v1/projects/4b53797444f846c4"},
            "size": 1187736145,
            "tags": [{"$ref": "/api/v1/tags/dbcf5b98bca54a16"}],
            "visibility": "public",
        },
        {
            "$uri": "/api/v1/samples/03242c0ab87048e1",
            "created_at": "2017-06-23T23:59:23.228640+00:00",
            "filename": "SRR4408292.fastq",
            "metadata": {"$ref": "/api/v1/metadata/3e7119ee74954abd"},
            "owner": {"$ref": "/api/v1/users/7189f36afe3640ac"},
            "primary_classification": {"$ref": "/api/v1/classifications/e0422602de41479f"},
            "project": {"$ref": "/api/v1/projects/4b53797444f846c4"},
            "size": 2267482973,
            "tags": [],
            "visibility": "public",
        },
    ],
    "GET::api/v1_experimental/functional_profiles\\?.*where=%7B%22sample%22%3A\\+%7B%22%24in%22%3A\\+%5B%2237e5151e7bcb4f87%22%5D%7D.*": [
        {
            "$uri": "/api/v1_experimental/functional_profiles/eec4ac90d9104d1e",
            "complete": True,
            "created_at": "2022-05-25T17:27:30.622286-07:00",
            "error_msg": "",
            "job": {"$ref": "/api/v1_experimental/jobs/59e7904ea8ed4202"},
            "sample": {"$ref": "/api/v1_experimental/samples/37e5151e7bcb4f87"},
            "success": True,
        },
    ],
    "GET::api/v1_experimental/functional_profiles\\?.*where=%7B%22sample%22%3A\\+%7B%22%24in%22%3A\\+%5B%2266c1531cb0b244f6%22%5D%7D.*": [
        {
            "$uri": "/api/v1_experimental/functional_profiles/bde18eb9407d4c2f",
            "complete": True,
            "created_at": "2022-05-25T17:27:30.622286-07:00",
            "error_msg": "",
            "job": {"$ref": "/api/v1_experimental/jobs/59e7904ea8ed4202"},
            "sample": {"$ref": "/api/v1_experimental/samples/66c1531cb0b244f6"},
            "success": True,
        },
    ],
    "GET::api/v1_experimental/functional_profiles\\?.*where=%7B%22sample%22%3A\\+%7B%22%24in%22%3A\\+%5B%22543c9c046e3e4e09%22%5D%7D.*": [
        {
            "$uri": "/api/v1_experimental/functional_profiles/31ddae978aff475f",
            "complete": True,
            "created_at": "2022-05-25T17:27:30.622286-07:00",
            "error_msg": "",
            "job": {"$ref": "/api/v1_experimental/jobs/59e7904ea8ed4202"},
            "sample": {"$ref": "/api/v1_experimental/samples/543c9c046e3e4e09"},
            "success": True,
        },
    ],
}

for functional_uuid in {"31ddae978aff475f", "bde18eb9407d4c2f", "eec4ac90d9104d1e"}:
    raw_results = json.load(
        open(
            f"tests/data/api/v1_experimental/functional_profiles/{functional_uuid}/results/index.json"
        )
    )
    for annotation in FunctionalAnnotations:
        for metric in FunctionalAnnotationsMetric.metrics_for_annotation(annotation):
            API_DATA[
                f"GET::api/v1_experimental/functional_profiles/{functional_uuid}/filtered_results\\?.*functional_group=%22{annotation}%22\\&metric=%22{metric}%22\\&taxa_stratified=true"
            ] = filtered_raw_results(raw_results, annotation, metric, True)
            API_DATA[
                f"GET::api/v1_experimental/functional_profiles/{functional_uuid}/filtered_results\\?.*functional_group=%22{annotation}%22\\&metric=%22{metric}%22\\&taxa_stratified=false"
            ] = filtered_raw_results(raw_results, annotation, metric, False)

SCHEMA_ROUTES = {}
API_DATA_DIR = os.path.join("tests", "data", "api")

for api_version in os.listdir(API_DATA_DIR):
    api_root = os.path.join(API_DATA_DIR, api_version)

    for resource_name in os.listdir(api_root):
        resource_root = os.path.join(api_root, resource_name)

        if resource_name == "schema":
            for filename in os.listdir(resource_root):
                if not filename.endswith(".json"):
                    continue

                if filename == "index.json":
                    resource_uri = f"GET::api/{api_version}/schema$"
                elif filename == "index_all.json":
                    resource_uri = f"GET::api/{api_version}/schema\\?expand=all"
                else:
                    resource_name = filename.replace(".json", "")
                    resource_uri = f"GET::api/{api_version}/{resource_name}/schema"

                SCHEMA_ROUTES[resource_uri] = json.load(open(os.path.join(resource_root, filename)))
        else:
            for dirpath, _, filenames in os.walk(resource_root):
                for filename in filenames:
                    if filename not in {"index.json", "index.json.gz"}:
                        continue

                    filepath = os.path.join(dirpath, filename)
                    if filepath.endswith(".json.gz"):
                        resource = json.load(gzip.open(filepath))
                    else:
                        resource = json.load(open(filepath))

                    # Remove tests/data/ from path and normalize path separator to "/" for the API
                    # route.
                    api_route = "/".join(dirpath.split(os.sep)[2:])
                    resource_uri = f"GET::{api_route}"
                    API_DATA[resource_uri] = resource

                    # If we're processing a resource's instance listing
                    # (e.g. api/v1/classifications), auto-mock each instance too
                    # (e.g. api/v1/classifications/<uuid>)
                    if dirpath == resource_root and isinstance(resource, list):
                        for instance in resource:
                            instance_uri = f"GET::{instance['$uri'].lstrip('/')}"
                            API_DATA[instance_uri] = instance


API_DATA.update(SCHEMA_ROUTES)


@pytest.yield_fixture(scope="function")
def api_data():
    with mock_requests(API_DATA) as rsps:
        yield rsps


@pytest.yield_fixture(scope="function")
def raw_api_data():
    yield copy.deepcopy(API_DATA)


@pytest.yield_fixture(scope="function")
def upload_mocks():
    def upload_callback(request):
        # Get and read the streaming iterator so it's empty
        if hasattr(request.body, "fields"):
            streaming_iterator = request.body.fields["file"][1]
            streaming_iterator.read()
        return (201, {"location": "on-aws"}, "")

    json_data = {
        "POST:multipart/form-data:fake_aws_callback": upload_callback,
        "GET::api/v1/samples/7428cca4a3a04a8e": {
            "$uri": "/api/v1/samples/7428cca4a3a04a8e",
            "created_at": "2015-09-25T17:27:19.596555-07:00",
            "filename": "SRR2352185.fastq.gz",
            "metadata": {"$ref": "/api/v1/metadata/a7fc7e430e704e2e"},
            "owner": {"$ref": "/api/v1/users/4ada56103d9a48b8"},
            "primary_classification": {"$ref": "/api/v1/classifications/464a7ebcf9f84050"},
            "project": None,
            "size": 181687821,
            "tags": [
                {"$ref": "/api/v1/tags/42997b7a62634985"},
                {"$ref": "/api/v1/tags/fb8e3b693c874f9e"},
                {"$ref": "/api/v1/tags/ff4e81909a4348d9"},
            ],
            "visibility": "private",
        },
        "POST::api/v1/samples/confirm_upload": "",
        "POST::api/v1/samples/init_multipart_upload": {
            "callback_url": "/api/import_file_from_s3",
            "file_id": "abcdef0987654321",
            "paired_end_file_id": "abcdefg0987654321",
            "s3_bucket": "onecodex-multipart-uploads-encrypted",
            "upload_aws_access_key_id": "aws_key",
            "upload_aws_secret_access_key": "aws_secret_key",
        },
        "POST::api/import_file_from_s3": "",
        "GET::api/v1/projects?.*where=%7B%22name.*": "",
        "GET::api/v1/projects?.*where=%7B%22project_name.*": [
            {
                "$uri": "/api/v1/projects/472fc57510e24150",
                "description": None,
                "name": "Test",
                "owner": {"$ref": "/api/v1/users/9923090af03c46ce"},
                "permissions": [
                    "can_see_files",
                    "can_incur_charges",
                    "can_download_files",
                    "can_edit_metadata",
                    "can_add_files",
                    "can_administer",
                ],
                "project_name": "testproj",
                "public": False,
            }
        ],
    }
    json_data.update(SCHEMA_ROUTES)
    with mock_requests(json_data):
        yield


# API FIXTURES
@pytest.yield_fixture(scope="session")
def ocx_schemas():
    with mock_requests(SCHEMA_ROUTES):
        yield


@pytest.fixture(scope="function")
def ocx():
    """Instantiated API client"""
    with mock_requests(SCHEMA_ROUTES):
        return Api(
            api_key="1eab4217d30d42849dbde0cd1bb94e39",
            base_url="http://localhost:3000",
            cache_schema=False,
        )


@pytest.fixture(scope="function")
def ocx_experimental():
    """Instantiated API client with experimental mode enabled"""
    with mock_requests(SCHEMA_ROUTES):
        return Api(
            api_key="1eab4217d30d42849dbde0cd1bb94e39",
            base_url="http://localhost:3000",
            cache_schema=False,
            experimental=True,
        )


@pytest.fixture(scope="function")
def custom_mock_requests():
    yield mock_requests


@pytest.fixture(scope="function")
def runner():
    runner = CliRunner(
        env={"ONE_CODEX_API_BASE": "http://localhost:3000", "ONE_CODEX_NO_TELEMETRY": "True"}
    )
    return runner


# SampleCollection fixtures
@pytest.fixture
def samples(ocx, api_data) -> SampleCollection:
    # Project with 3 samples
    return ocx.Samples.where(project="4b53797444f846c4")


# CLI / FILE SYSTEM FIXTURE
@pytest.fixture(scope="function")
def mocked_creds_path(monkeypatch, tmpdir):
    # TODO: tmpdir is actually a LocalPath object
    # from py.path, and we coerce it into a string
    # for compatibility with the existing library code
    # *but* we should perhaps *not* do that for
    # better cross-platform compatibility. Investigate
    # and update as needed.
    def mockreturn(path):
        return os.path.join(str(tmpdir), ".onecodex")

    monkeypatch.setattr(os.path, "expanduser", mockreturn)


@pytest.fixture(scope="function")
def mocked_creds_file(mocked_creds_path):
    with open(os.path.expanduser("~/.onecodex"), mode="w") as f:
        now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

        f.write(
            json.dumps(
                {
                    "email": "test@onecodex.com",
                    "api_key": "123yuixha87yd87q3123uiqhsd8q2738",
                    "saved_at": now,
                    "updated_at": now,
                }
            )
        )


# Upload mock helpers
FASTQ_SEQUENCE = """
@cluster_2:UMI_ATTCCG
TTTCCGGGGCACATAATCTTCAGCCGGGCGC
+
9C;=;=<9@4868>9:67AA<9>65<=>591"""


@contextmanager
def path_for_filename(tmp_path, runner, filename):
    path = os.path.join(str(tmp_path), filename)
    parent_dir = os.path.dirname(path)
    with runner.isolated_filesystem():
        if not os.path.exists(parent_dir):
            os.makedirs(parent_dir)
        yield path


@pytest.fixture
def generate_fastq(tmp_path, runner):
    def fn(filename):
        with path_for_filename(tmp_path, runner, filename) as path:
            with open(path, "w") as fout:
                fout.write(FASTQ_SEQUENCE)
            return path

    yield fn


@pytest.fixture
def generate_fastq_gz(tmp_path, runner):
    def fn(filename):
        with path_for_filename(tmp_path, runner, filename) as path:
            with gzip.open(path, "w") as fout:
                fout.write(FASTQ_SEQUENCE.encode("utf-8"))
            return path

    yield fn
