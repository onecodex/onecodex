import json
import os
import os.path

import mock
import pytest
from click.testing import CliRunner
from testfixtures import Replace

from onecodex import Cli
from tests.conftest import API_DATA

DATE_FORMAT = "%Y-%m-%d %H:%M"


@pytest.fixture
def mock_file_upload():
    with mock.patch("onecodex.lib.upload._upload_sequence_fileobj") as m:
        yield m


@pytest.fixture
def mock_sample_get():
    with mock.patch("onecodex.models.sample.Samples.get") as m:
        yield m


# Click 8.2.0 changed exit code from 0 to 2 if no arguments are supplied
@pytest.mark.parametrize("args,exit_codes", [(None, {0, 2}), ("-h", {0}), ("--help", {0})])
def test_cli_help(runner, args, exit_codes):
    command = [args] if args is not None else None
    result = runner.invoke(Cli, command)
    assert result.exit_code in exit_codes
    assert "One Codex v1 API command line interface" in result.output


def test_version(runner):
    result = runner.invoke(Cli, ["--version"])
    assert result.exit_code == 0
    assert "onecodex, version" in result.output


# Test CLI without base override
def test_cli_wo_override(api_data, monkeypatch):
    monkeypatch.delattr("requests.sessions.Session.request")
    runner = CliRunner()
    result = runner.invoke(Cli, ["-v", "--help"])
    assert result.exit_code == 0
    assert "One Codex v1 API command line interface" in result.output


# Analyses
def test_analysis_help(runner, api_data, mocked_creds_path):
    result = runner.invoke(Cli, ["analyses", "--help"])
    analysis_desc = "Retrieve performed analyses"
    assert result.exit_code == 0
    assert analysis_desc in result.output


def test_analyses(runner, api_data, mocked_creds_file):
    r0 = runner.invoke(Cli, ["analyses"])
    r1 = runner.invoke(Cli, ["analyses", "593601a797914cbf"])
    assert r0.exit_code == 0
    assert r1.exit_code == 0
    assert API_DATA["GET::api/v1/analyses/593601a797914cbf"]["$uri"] in r0.output
    assert API_DATA["GET::api/v1/analyses/593601a797914cbf"]["$uri"] in r1.output


def _analysis_payload(analysis_id, **overrides):
    payload = {
        "$uri": f"/api/v1/analyses/{analysis_id}",
        "analysis_type": "classification",
        "created_at": "2015-09-25T17:27:30.622286-07:00",
        "complete": True,
        "success": True,
        "error_msg": None,
        "job": {"$ref": "/api/v1/jobs/e4b1ab37ff554c53"},
        "sample": {"$ref": "/api/v1/samples/7428cca4a3a04a8e"},
        "cost": None,
        "dependencies": [],
        "draft": False,
        "job_args": {},
    }
    payload.update(overrides)
    return payload


def test_analyses_await(runner, custom_mock_requests, mocked_creds_file, monkeypatch):
    analysis_id = "593601a797914cbf"
    bodies = [
        _analysis_payload(analysis_id, complete=False, success=False),
        _analysis_payload(analysis_id, complete=False, success=False),
        _analysis_payload(analysis_id),
    ]

    def get_callback(request):
        body = bodies.pop(0)
        return (200, {"Content-Type": "application/json"}, json.dumps(body))

    with custom_mock_requests({f"GET::api/v1/analyses/{analysis_id}": get_callback}):
        with mock.patch("onecodex.models.analysis.time.sleep"):
            result = runner.invoke(Cli, ["await", analysis_id, "--initial-interval", "5"])

    assert result.exit_code == 0, result.output
    assert f"Analysis {analysis_id} complete" in result.output
    assert "success=True" in result.output


def test_analyses_logs(runner, custom_mock_requests, mocked_creds_file):
    analysis_id = "593601a797914cbf"

    captured = {}

    def get_callback(request):
        return (
            200,
            {"Content-Type": "application/json"},
            json.dumps(_analysis_payload(analysis_id)),
        )

    def logs_callback(request):
        captured["url"] = request.url
        return (200, {"Content-Type": "text/plain"}, "log line A\nlog line B\n")

    with custom_mock_requests(
        {
            f"GET::api/v1/analyses/{analysis_id}": get_callback,
            f"GET:text/plain:api/v1/analyses/{analysis_id}/logs": logs_callback,
        }
    ):
        result = runner.invoke(Cli, ["analyses", "logs", analysis_id, "--tail", "2"])

    assert result.exit_code == 0, result.output
    assert "log line A" in result.output
    assert "log line B" in result.output
    assert "tail=2" in captured["url"]


def test_analyses_logs_404(runner, custom_mock_requests, mocked_creds_file):
    analysis_id = "593601a797914cbf"

    def get_callback(request):
        return (
            200,
            {"Content-Type": "application/json"},
            json.dumps(_analysis_payload(analysis_id)),
        )

    def logs_callback(request):
        return (404, {"Content-Type": "text/plain"}, "Not Found")

    with custom_mock_requests(
        {
            f"GET::api/v1/analyses/{analysis_id}": get_callback,
            f"GET:text/plain:api/v1/analyses/{analysis_id}/logs": logs_callback,
        }
    ):
        result = runner.invoke(Cli, ["analyses", "logs", analysis_id])

    assert result.exit_code != 0
    assert "Logs not found" in result.output
    assert "Traceback" not in result.output


FILE_DETAILS = {
    "files": [
        {
            "filename": "report.tsv",
            "filepath": "results/report.tsv",
            "size": 11,
            "url": "http://localhost:3000/files/report.tsv",
        },
        {
            "filename": "summary.txt",
            "filepath": "summary.txt",
            "size": 7,
            "url": "http://localhost:3000/files/summary.txt",
        },
    ]
}


@pytest.fixture
def mock_analysis_files(custom_mock_requests):
    from onecodex.models.analysis import _AnalysesBase

    _AnalysesBase.get_files.cache_clear()

    analysis_id = "593601a797914cbf"

    def make_mocks(file_details=FILE_DETAILS):
        return custom_mock_requests(
            {
                f"GET::api/v1/analyses/{analysis_id}": _analysis_payload(analysis_id),
                f"GET::api/v1/analyses/{analysis_id}/file_details": file_details,
                "GET:text/plain:files/report.tsv": lambda request: (
                    200,
                    {"Content-Type": "text/plain"},
                    "report data",
                ),
                "GET:text/plain:files/summary.txt": lambda request: (
                    200,
                    {"Content-Type": "text/plain"},
                    "summary",
                ),
            }
        )

    return analysis_id, make_mocks


def test_analyses_files(runner, mock_analysis_files, mocked_creds_file):
    analysis_id, make_mocks = mock_analysis_files
    with make_mocks():
        result = runner.invoke(Cli, ["analyses", "files", analysis_id])

    assert result.exit_code == 0, result.output
    assert [line.split() for line in result.output.splitlines()] == [
        ["Filepath", "Size"],
        ["-" * 66, "-" * 9],
        ["results/report.tsv", "11", "B"],
        ["summary.txt", "7", "B"],
    ]


def test_analyses_files_long_filepath_and_size(runner, mock_analysis_files, mocked_creds_file):
    analysis_id, make_mocks = mock_analysis_files
    file_details = {
        "files": [
            {
                "filename": "report.tsv",
                "filepath": "results/" + "nested/" * 10 + "report.tsv",
                "size": 12300000,
                "url": "http://localhost:3000/files/report.tsv",
            }
        ]
    }
    with make_mocks(file_details=file_details):
        result = runner.invoke(Cli, ["analyses", "files", analysis_id])

    assert result.exit_code == 0, result.output
    filepath, size, unit = result.output.splitlines()[2].split()
    assert filepath == "results/" + "nested/" * 7 + "nested..."
    assert len(filepath) == 66
    assert (size, unit) == ("12.3", "MB")


def test_analyses_files_json(runner, mock_analysis_files, mocked_creds_file):
    analysis_id, make_mocks = mock_analysis_files
    with make_mocks():
        result = runner.invoke(Cli, ["analyses", "files", analysis_id, "--json"])

    assert result.exit_code == 0, result.output
    assert json.loads(result.output) == FILE_DETAILS["files"]


def test_analyses_files_empty(runner, mock_analysis_files, mocked_creds_file):
    analysis_id, make_mocks = mock_analysis_files
    with make_mocks(file_details={"files": []}):
        result = runner.invoke(Cli, ["analyses", "files", analysis_id])

    assert result.exit_code == 0, result.output
    assert "no output files" in result.output


def test_analyses_download(runner, mock_analysis_files, mocked_creds_file):
    analysis_id, make_mocks = mock_analysis_files
    with runner.isolated_filesystem():
        with make_mocks():
            result = runner.invoke(Cli, ["analyses", "download", analysis_id, "-o", "out"])

        assert result.exit_code == 0, result.output
        with open(os.path.join("out", "results", "report.tsv")) as f:
            assert f.read() == "report data"
        with open(os.path.join("out", "summary.txt")) as f:
            assert f.read() == "summary"


def test_analyses_download_single_file(runner, mock_analysis_files, mocked_creds_file):
    analysis_id, make_mocks = mock_analysis_files
    with runner.isolated_filesystem():
        with make_mocks():
            result = runner.invoke(
                Cli, ["analyses", "download", analysis_id, "-f", "results/report.tsv"]
            )

        assert result.exit_code == 0, result.output
        assert os.path.exists(os.path.join("results", "report.tsv"))
        assert not os.path.exists("summary.txt")


def test_analyses_download_nonexistent_file(runner, mock_analysis_files, mocked_creds_file):
    analysis_id, make_mocks = mock_analysis_files
    with runner.isolated_filesystem():
        with make_mocks():
            result = runner.invoke(Cli, ["analyses", "download", analysis_id, "-f", "nope.tsv"])

        assert result.exit_code != 0
        assert "has no output file(s): nope.tsv" in result.output


def test_analyses_download_refuses_path_traversal(runner, mock_analysis_files, mocked_creds_file):
    analysis_id, make_mocks = mock_analysis_files
    file_details = {
        "files": [
            {
                "filename": "report.tsv",
                "filepath": "../report.tsv",
                "size": 11,
                "url": "http://localhost:3000/files/report.tsv",
            }
        ]
    }
    with runner.isolated_filesystem():
        with make_mocks(file_details=file_details):
            result = runner.invoke(Cli, ["analyses", "download", analysis_id, "-o", "out"])

        assert result.exit_code != 0
        assert "Refusing to write '../report.tsv'" in result.output
        assert not os.path.exists(os.path.join("..", "report.tsv"))


@pytest.mark.parametrize(
    "args",
    [
        ["analyses", "--help"],
        ["analyses", "files", "--help"],
        ["analyses", "download", "--help"],
        ["analyses", "logs", "--help"],
        ["analyses", "await", "--help"],
    ],
)
def test_analyses_help_does_not_require_login(runner, mocked_creds_path, args):
    result = runner.invoke(Cli, args)

    assert result.exit_code == 0, result.output
    assert "requires authentication" not in result.output


def test_analyses_download_existing_file_aborts(runner, mock_analysis_files, mocked_creds_file):
    analysis_id, make_mocks = mock_analysis_files
    with runner.isolated_filesystem():
        os.makedirs(os.path.join("out", "results"))
        with open(os.path.join("out", "results", "report.tsv"), "w") as f:
            f.write("mine")

        with make_mocks():
            result = runner.invoke(Cli, ["analyses", "download", analysis_id, "-o", "out"])

        assert result.exit_code != 0
        assert "already exists. Will not overwrite." in result.output
        with open(os.path.join("out", "results", "report.tsv")) as f:
            assert f.read() == "mine"
        assert not os.path.exists(os.path.join("out", "summary.txt"))


def test_analyses_download_existing_outdir(runner, mock_analysis_files, mocked_creds_file):
    analysis_id, make_mocks = mock_analysis_files
    with runner.isolated_filesystem():
        os.makedirs("out")
        with make_mocks():
            result = runner.invoke(Cli, ["analyses", "download", analysis_id, "-o", "out"])

        assert result.exit_code == 0, result.output
        assert os.path.exists(os.path.join("out", "results", "report.tsv"))


def test_analyses_download_outdir_is_a_file(runner, mock_analysis_files, mocked_creds_file):
    analysis_id, make_mocks = mock_analysis_files
    with runner.isolated_filesystem():
        with open("out", "w") as f:
            f.write("not a directory")

        with make_mocks():
            result = runner.invoke(Cli, ["analyses", "download", analysis_id, "-o", "out"])

        assert result.exit_code != 0
        assert "is a file" in result.output


def _file_details(*filepaths):
    return {
        "files": [
            {
                "filename": os.path.basename(filepath),
                "filepath": filepath,
                "size": 11,
                "url": "http://localhost:3000/files/report.tsv",
            }
            for filepath in filepaths
        ]
    }


def test_analyses_download_file_and_directory_collide(
    runner, mock_analysis_files, mocked_creds_file
):
    analysis_id, make_mocks = mock_analysis_files
    with runner.isolated_filesystem():
        with make_mocks(file_details=_file_details("a", "a/b")):
            result = runner.invoke(Cli, ["analyses", "download", analysis_id, "-o", "out"])

        assert result.exit_code != 0
        assert "lists 'a' as an output file and as a directory containing 'a/b'" in result.output
        assert "Traceback" not in result.output
        assert not os.path.exists("out")


def test_analyses_download_directory_blocked_by_existing_file(
    runner, mock_analysis_files, mocked_creds_file
):
    analysis_id, make_mocks = mock_analysis_files
    with runner.isolated_filesystem():
        os.makedirs("out")
        with open(os.path.join("out", "a"), "w") as f:
            f.write("mine")

        with make_mocks(file_details=_file_details("a/b/c")):
            result = runner.invoke(Cli, ["analyses", "download", analysis_id, "-o", "out"])

        assert result.exit_code != 0
        assert "already exists and is not a directory" in result.output
        assert "Traceback" not in result.output
        with open(os.path.join("out", "a")) as f:
            assert f.read() == "mine"


def test_analyses_download_duplicate_filepaths(runner, mock_analysis_files, mocked_creds_file):
    analysis_id, make_mocks = mock_analysis_files
    with runner.isolated_filesystem():
        with make_mocks(file_details=_file_details("a/b", "/a/b")):
            result = runner.invoke(Cli, ["analyses", "download", analysis_id, "-o", "out"])

        assert result.exit_code != 0
        assert "which both resolve to" in result.output
        assert not os.path.exists("out")


@pytest.mark.parametrize("filepath", ["", "/", "."])
def test_analyses_download_degenerate_filepath(
    runner, mock_analysis_files, mocked_creds_file, filepath
):
    analysis_id, make_mocks = mock_analysis_files
    with runner.isolated_filesystem():
        with make_mocks(file_details=_file_details(filepath)):
            result = runner.invoke(Cli, ["analyses", "download", analysis_id, "-o", "out"])

        assert result.exit_code != 0
        assert "Refusing to write" in result.output


# Classifications
def test_classification_instance(runner, api_data, mocked_creds_file):
    result = runner.invoke(Cli, ["classifications", "593601a797914cbf"])
    assert result.exit_code == 0
    assert API_DATA["GET::api/v1/classifications/593601a797914cbf"]["$uri"] in result.output


def test_classifications_table(runner, api_data, mocked_creds_file, monkeypatch):
    result = runner.invoke(Cli, ["classifications", "45a573fb7833449a", "--results"])
    assert result.exit_code == 0
    assert "Staphylococcus" in result.output


# Documents
def test_documents_table(runner, api_data, mocked_creds_file):
    result = runner.invoke(Cli, ["documents", "list"])
    assert len(result.stdout.split("\n")) == 20
    assert "OneCodexTakeHome" in result.stdout

    result = runner.invoke(Cli, ["documents", "list", "--json"])
    assert len(result.stdout.split("\n")) == 258  # shorter because we ignore extra fields
    assert "OneCodexTakeHome" in result.stdout


# Assets
def test_assets_table(runner, api_data, mocked_creds_file):
    result = runner.invoke(Cli, ["assets", "list"])
    assert result.exit_code == 0
    assert "ID" in result.stdout
    assert "Status" in result.stdout
    assert "057268a2ba9f6d7a" in result.stdout
    assert "reference_genome.fa.gz" in result.stdout
    assert "available" in result.stdout
    assert "importing" in result.stdout

    result = runner.invoke(Cli, ["assets", "list", "--json"])
    assert result.exit_code == 0
    assert "057268a2ba9f6d7a" in result.stdout
    assert "reference_genome.fa.gz" in result.stdout


def test_assets_list_empty(runner, custom_mock_requests, mocked_creds_file):
    with custom_mock_requests({"GET::api/v1/assets": []}):
        result = runner.invoke(Cli, ["assets", "list"])
        assert result.exit_code == 0
        assert "haven't uploaded any assets" in result.stdout


# Panels
def test_panel_instances(runner, api_data, mocked_creds_file):
    result = runner.invoke(Cli, ["panels"])
    assert result.exit_code == 0


# Workflows
def test_workflow_instances(runner, api_data, mocked_creds_file):
    result = runner.invoke(Cli, ["workflows"])
    assert result.exit_code == 0


# Jobs create/update
def test_jobs_create_cli(runner, api_data, custom_mock_requests, mocked_creds_file, tmp_path):
    new_job_id = "0123456789abcdef"
    captured = {}

    def create_callback(request):
        captured["body"] = json.loads(request.body)
        return (
            200,
            {"Content-Type": "application/json"},
            json.dumps(
                {
                    "$uri": f"/api/v1/jobs/{new_job_id}",
                    "created_at": "2026-01-01T00:00:00+00:00",
                    "name": "my-job",
                    "analysis_type": "custom",
                    "public": False,
                    "draft": False,
                    "job_args_schema": {},
                }
            ),
        )

    script = tmp_path / "run.sh"
    script.write_text("echo hello\n")
    schema = tmp_path / "schema.json"
    schema.write_text(json.dumps([{"fields": [{"name": "min_quality", "type": "integer"}]}]))

    asset_id = "057268a2ba9f6d7a"
    parent_job_id = "cc1d331e1ee54bac"

    with custom_mock_requests({"POST::api/v1/jobs": create_callback}):
        result = runner.invoke(
            Cli,
            [
                "jobs",
                "create",
                *("--name", "my-job"),
                *("--script", str(script)),
                *("--image-uri", "docker.io/library/python:3.12"),
                *("--cpu", "1"),
                *("--ram-gb", "1"),
                *("--storage-gb", "1"),
                *("--asset-id", asset_id),
                *("-d", f"{parent_job_id}=outdir"),
                *("--arguments-schema", str(schema)),
            ],
        )

    assert result.exit_code == 0, result.output
    assert f"Created job {new_job_id}" in result.output
    assert captured["body"] == {
        "name": "my-job",
        "script": "echo hello\n",
        "image_uri": "docker.io/library/python:3.12",
        "cpu": 1.0,
        "ram_gb": 1.0,
        "storage_gb": 1.0,
        "assets": [{"$ref": f"/api/v1/assets/{asset_id}"}],
        "dependencies": [
            {"job": {"$ref": f"/api/v1/jobs/{parent_job_id}"}, "output_dir": "outdir"}
        ],
        "arguments_schema": [{"fields": [{"name": "min_quality", "type": "integer"}]}],
    }


def _jobs_write_mock(captured, job_id="0123456789abcdef"):
    def callback(request):
        captured["body"] = json.loads(request.body)
        return (
            200,
            {"Content-Type": "application/json"},
            json.dumps(
                {
                    "$uri": f"/api/v1/jobs/{job_id}",
                    "created_at": "2026-01-01T00:00:00+00:00",
                    "name": "my-nf-job",
                    "analysis_type": "workflow",
                    "public": False,
                    "draft": True,
                    "job_args_schema": {},
                }
            ),
        )

    return callback


def test_jobs_create_nextflow_cli(
    runner, api_data, custom_mock_requests, mocked_creds_file, tmp_path
):
    new_job_id = "0123456789abcdef"
    captured = {}

    script = tmp_path / "run_nextflow.sh"
    script.write_text("nextflow run main.nf\n")

    with custom_mock_requests({"POST::api/v1/jobs": _jobs_write_mock(captured, new_job_id)}):
        result = runner.invoke(
            Cli,
            [
                "jobs",
                "create",
                *("--name", "my-nf-job"),
                *("--script", str(script)),
                *("--job-type", "nextflow"),
                *("--nextflow-version", "25.10.0"),
                *("--repository-url", "https://github.com/org/repo"),
                *("--repository-tag", "v0.1.0"),
            ],
        )

    assert result.exit_code == 0, result.output
    assert captured["body"] == {
        "name": "my-nf-job",
        "script": "nextflow run main.nf\n",
        "nextflow_version": "25.10.0",
        "job_type": "nextflow",
        "repository": {"url": "https://github.com/org/repo", "tag": "v0.1.0"},
    }


def test_jobs_update_nextflow_version_cli(
    runner, api_data, custom_mock_requests, mocked_creds_file
):
    job_id = "cc1d331e1ee54bac"
    captured = {}

    with custom_mock_requests({f"PATCH::api/v1/jobs/{job_id}": _jobs_write_mock(captured, job_id)}):
        result = runner.invoke(Cli, ["jobs", "update", job_id, *("--nextflow-version", "24.10.5")])

    assert result.exit_code == 0, result.output
    assert captured["body"] == {"nextflow_version": "24.10.5"}


@pytest.mark.parametrize(
    "image_args,expected_error",
    [
        ([], "image-uri"),
        (
            ["--job-type", "nextflow", "--image-uri", "docker.io/library/python:3.12"],
            "--image-uri is not supported for Nextflow jobs",
        ),
        (
            ["--image-uri", "docker.io/library/python:3.12", "--nextflow-version", "25.10.0"],
            "--nextflow-version requires --job-type nextflow",
        ),
    ],
)
def test_jobs_create_invalid_image_options(
    runner, mocked_creds_file, tmp_path, image_args, expected_error
):
    script = tmp_path / "run.sh"
    script.write_text("echo hi\n")
    result = runner.invoke(
        Cli, ["jobs", "create", *("--name", "x"), *("--script", str(script)), *image_args]
    )
    assert result.exit_code != 0
    assert expected_error in result.output


def test_jobs_update_cli(runner, api_data, custom_mock_requests, mocked_creds_file):
    job_id = "cc1d331e1ee54bac"
    captured = {}

    def patch_callback(request):
        captured["body"] = json.loads(request.body)
        return (
            200,
            {"Content-Type": "application/json"},
            json.dumps(
                {
                    "$uri": f"/api/v1/jobs/{job_id}",
                    "created_at": "2016-05-05T17:27:02.116480+00:00",
                    "name": "Renamed",
                    "analysis_type": "classification",
                    "public": True,
                    "draft": False,
                    "job_args_schema": {},
                }
            ),
        )

    with custom_mock_requests({f"PATCH::api/v1/jobs/{job_id}": patch_callback}):
        result = runner.invoke(Cli, ["jobs", "update", job_id, *("--name", "Renamed")])

    assert result.exit_code == 0, result.output
    assert captured["body"] == {"name": "Renamed"}
    assert f"Updated job {job_id}" in result.output


def _jobs_run_mock(captured, analysis_id):
    def callback(request):
        captured["body"] = json.loads(request.body)
        return (
            200,
            {"Content-Type": "application/json"},
            json.dumps({"$ref": f"/api/v1/analyses/{analysis_id}"}),
        )

    return callback


def test_jobs_run_args_json(runner, api_data, custom_mock_requests, mocked_creds_file):
    job_id = "47c4fe23588640a9"
    sample_id = "7428cca4a3a04a8e"
    analysis_id = "593601a797914cbf"
    captured = {}

    with custom_mock_requests(
        {f"POST::api/v1/jobs/{job_id}/run": _jobs_run_mock(captured, analysis_id)}
    ):
        result = runner.invoke(
            Cli,
            [
                "jobs",
                "run",
                job_id,
                sample_id,
                *("--args-json", '{"min_quality": 30, "trim": true, "adapter": "AGATC"}'),
            ],
        )

    assert result.exit_code == 0, result.output
    assert captured["body"]["job_args"] == {
        "min_quality": 30,
        "trim": True,
        "adapter": "AGATC",
    }


def test_jobs_run_args_json_mutually_exclusive(runner, mocked_creds_file):
    result = runner.invoke(
        Cli,
        [
            "jobs",
            "run",
            "47c4fe23588640a9",
            "7428cca4a3a04a8e",
            *("-a", "x=1"),
            *("--args-json", "{}"),
        ],
    )
    assert result.exit_code != 0
    assert "Cannot combine" in result.output


def test_jobs_run_args_json_invalid_json(runner, mocked_creds_file):
    result = runner.invoke(
        Cli,
        ["jobs", "run", "47c4fe23588640a9", "7428cca4a3a04a8e", *("--args-json", "not json")],
    )
    assert result.exit_code != 0
    assert "Invalid JSON" in result.output


def test_jobs_run_args_json_not_object(runner, mocked_creds_file):
    result = runner.invoke(
        Cli,
        ["jobs", "run", "47c4fe23588640a9", "7428cca4a3a04a8e", *("--args-json", "[1, 2, 3]")],
    )
    assert result.exit_code != 0
    assert "JSON object" in result.output


def test_jobs_update_no_fields(runner, api_data, mocked_creds_file):
    job_id = "cc1d331e1ee54bac"
    result = runner.invoke(Cli, ["jobs", "update", job_id])
    assert result.exit_code != 0
    assert "No fields provided" in result.output


def test_jobs_update_rejects_job_type(runner, api_data, mocked_creds_file):
    job_id = "cc1d331e1ee54bac"
    result = runner.invoke(Cli, ["jobs", "update", job_id, *("--job-type", "shell_script")])
    assert result.exit_code != 0
    assert "job-type" in result.output


# Samples
def test_samples(runner, api_data, mocked_creds_file):
    r0 = runner.invoke(Cli, ["samples"])
    r1 = runner.invoke(Cli, ["samples", "7428cca4a3a04a8e"])
    assert r0.exit_code == 0
    assert r1.exit_code == 0
    assert API_DATA["GET::api/v1/samples/7428cca4a3a04a8e"]["$uri"] in r0.output
    assert API_DATA["GET::api/v1/samples/7428cca4a3a04a8e"]["$uri"] in r1.output


# Login tests
def mock_fetch_api_key(username, password, server_url):
    return "123yuixha87yd87q3123uiqhsd8q2738"


def test_api_login(runner, mocked_creds_path):
    login_input = "user@example.com" + "\n" + "userpassword" + "\n"
    successful_login_msg = "Your ~/.onecodex credentials file was successfully created."
    with Replace("onecodex.auth.fetch_api_key_from_uname", mock_fetch_api_key):
        result = runner.invoke(Cli, ["login"], input=login_input)
        assert result.exit_code == 0
        assert successful_login_msg in result.output


@pytest.mark.parametrize(
    "email,success,code",
    [("incorrect+email@onecodex.com", False, 1), ("asmngs@onecodex.com", True, 0)],
)
def test_api_key_login(runner, api_data, mocked_creds_path, email, success, code):
    result = runner.invoke(Cli, ["--api-key", "0" * 32, "login"], input=email)
    assert result.exit_code == code
    assert os.path.exists(os.path.expanduser("~/.onecodex")) is success


def test_creds_file_exists(runner, mocked_creds_file):
    expected_message = "Credentials file already exists"

    result = runner.invoke(Cli, ["login"])
    assert result.exit_code == 0
    assert expected_message in result.output


def test_silent_login(runner, mocked_creds_file, api_data):
    result = runner.invoke(Cli, ["samples"])
    assert result.exit_code == 0


def test_creds_file_corrupted(runner, mocked_creds_file):
    path = os.path.expanduser("~/.onecodex")
    with open(path, mode="w") as f:
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
    result = runner.invoke(Cli, ["samples"], catch_exceptions=False)
    assert "requires authentication" in result.output

    # bearer token in environment
    with mock.patch.dict(
        os.environ, {"ONE_CODEX_BEARER_TOKEN": "00000000000000000000000000000000"}
    ):
        assert "ONE_CODEX_BEARER_TOKEN" in os.environ
        assert "ONE_CODEX_API_KEY" not in os.environ
        result = runner.invoke(Cli, ["samples"], catch_exceptions=False)
    assert "requires authentication" not in result.output

    # bearer token in environment
    with mock.patch.dict(os.environ, {"ONE_CODEX_API_KEY": "00000000000000000000000000000000"}):
        assert "ONE_CODEX_BEARER_TOKEN" not in os.environ
        assert "ONE_CODEX_API_KEY" in os.environ
        result = runner.invoke(Cli, ["samples"], catch_exceptions=False)
    assert "requires authentication" not in result.output


# Uploads
@pytest.mark.parametrize(
    "files,threads",
    [
        (["temp.fa"], False),
        (["temp.fa"], True),
        (["temp.fa", "temp2.fa"], False),
        (["temp.fa", "temp2.fa"], True),
    ],
)
def test_standard_uploads(
    runner, mocked_creds_path, caplog, generate_fastq, upload_mocks, files, threads
):
    """Test single and multi file uploads, with and without threads
    (but not files >5GB)
    """
    with (
        runner.isolated_filesystem(),
        mock.patch("onecodex.models.Projects.get", side_effect=lambda _: None),
        mock.patch("onecodex.lib.upload._s3_intermediate_upload") as s3_upload_mock,
    ):
        s3_upload_mock.return_value = {"sample_id": "7428cca4a3a04a8e"}
        args = [
            "--api-key",
            "01234567890123456789012345678901",
            "-v",
            "upload",
            "--project",
            "testproj",
            "--tag",
            "isolate",
            "--metadata",
            "totalige=1234",
        ]
        if not threads:
            args += ["--max-threads", "1"]
        for f in files:
            args.append(generate_fastq(f))

        result = runner.invoke(Cli, args, catch_exceptions=False, input="\n")
        assert result.exit_code == 0
        assert "7428cca4a3a04a8e" in caplog.text  # mocked file id


def test_empty_upload(runner, mocked_creds_path, upload_mocks):
    with runner.isolated_filesystem():
        f = "tmp.fa"
        f_out = open(f, mode="w")
        f_out.close()
        args = ["--api-key", "01234567890123456789012345678901", "upload", f]
        result = runner.invoke(Cli, args, catch_exceptions=False)
        assert result.exit_code != 0


@pytest.mark.parametrize(
    "files,n_samples_uploaded,n_files_uploaded,n_paired_files,n_multiline_groups",
    [
        # 1 files, 1 sample
        (["test.fq"], 1, 1, 0, 0),
        # 2 files, 1 sample
        (["test_R1.fq", "test_R2.fq"], 1, 2, 2, 0),
        (["test_1.fq", "test_2.fq"], 1, 2, 2, 0),
        # 2 files, 2 lines each, 1 sample
        (
            [
                "test_S1_L001_R1_001.fastq.gz",
                "test_S1_L001_R2_001.fastq.gz",
                "test_S1_L002_R1_001.fastq.gz",
                "test_S1_L002_R2_001.fastq.gz",
            ],
            1,
            2,
            4,
            1,
        ),
        (["dir_r1_test/test_R1.fq", "dir_r1_test/test_R2.fq"], 1, 2, 2, 0),
        (["test_1_1.fq", "test_1_2.fq"], 1, 2, 2, 0),
        # 3 files, 2 samples
        (["test_R1.fq", "test_R2.fq", "other.fq"], 2, 3, 2, 0),
    ],
)
def test_paired_and_multiline_files(
    runner,
    generate_fastq,
    mock_file_upload,
    mock_sample_get,
    mocked_creds_path,
    upload_mocks,
    files,
    n_samples_uploaded,
    n_files_uploaded,
    n_paired_files,
    n_multiline_groups,
):
    files = [generate_fastq(x) for x in files]

    args = ["--api-key", "01234567890123456789012345678901", "upload"] + files
    # check that 2 uploads are kicked off for the pair of files
    result = runner.invoke(Cli, args, catch_exceptions=False, input="\n\n\n")
    assert mock_file_upload.call_count == n_files_uploaded
    assert mock_sample_get.call_count == n_samples_uploaded
    assert result.exit_code == 0

    paired_files_prompt = "It appears there are {} paired files (of {} total)".format(
        n_paired_files, len(files)
    )
    if n_paired_files > 0:
        assert paired_files_prompt in result.output
    else:
        assert paired_files_prompt not in result.output

    multilane_prompt = "This data appears to have been split across multiple sequencing lanes.\nConcatenate lanes before upload?"
    if n_multiline_groups > 0:
        assert multilane_prompt in result.output
    else:
        assert multilane_prompt not in result.output


def test_paired_files_with_forward_and_reverse_args(
    runner, generate_fastq, mock_file_upload, mock_sample_get, mocked_creds_path, upload_mocks
):
    f1, f2 = generate_fastq("test_R1.fq"), generate_fastq("test_R2.fq")

    # Check with --forward and --reverse, should succeed
    args = [
        "--api-key",
        "01234567890123456789012345678901",
        "upload",
        "--forward",
        f1,
        "--reverse",
        f2,
    ]
    result = runner.invoke(Cli, args, input="Y")
    assert mock_file_upload.call_count == 2
    assert mock_sample_get.call_count == 1
    assert "It appears there are 2 paired files" not in result.output  # skips message
    assert result.exit_code == 0

    # Check with only --forward, should fail
    args = ["--api-key", "01234567890123456789012345678901", "upload", "--forward", f1]
    result = runner.invoke(Cli, args)
    assert "You must specify both forward and reverse files" in result.output
    assert result.exit_code != 0

    # Check with only --reverse, should fail
    args = ["--api-key", "01234567890123456789012345678901", "upload", "--reverse", f2]
    result = runner.invoke(Cli, args)
    assert "You must specify both forward and reverse files" in result.output
    assert result.exit_code != 0

    # Check with --forward, --reverse, and files, should fail
    args = [
        "--api-key",
        "01234567890123456789012345678901",
        "upload",
        "--forward",
        f1,
        "--reverse",
        f2,
        f2,
    ]
    result = runner.invoke(Cli, args)
    assert "You may not pass a FILES argument" in result.output
    assert result.exit_code != 0


@pytest.mark.parametrize(
    "files,n_samples_uploaded,n_files_uploaded,n_paired_files,n_ont_files",
    [
        # 1 files, 1 sample
        (["test.fq"], 1, 1, 0, 0),
        # 2 files, 1 sample
        (["test_1.fq", "test_2.fq"], 1, 2, 2, 0),
        # 4 files, 4 ONT parts, 1 sample
        (
            [
                "test_S1_0.fastq.gz",
                "test_S1_1.fastq.gz",
                "test_S1_2.fastq.gz",
                "test_S1_3.fastq.gz",
            ],
            1,
            1,
            0,
            4,
        ),
        # 3 files, 2 ONT parts, 2 samples
        (
            [
                "test_S1_0.fastq.gz",
                "test_S1_1.fastq.gz",
                "test_S2_0.fastq.gz",
            ],
            2,
            2,
            0,
            2,
        ),
        # 6 files, 2x3 ONT parts, 2 samples
        (
            [
                "test_S1_0.fastq.gz",
                "test_S2_0.fastq.gz",
                "test_S1_1.fastq.gz",
                "test_S2_1.fastq.gz",
                "test_S1_2.fastq.gz",
                "test_S2_2.fastq.gz",
            ],
            2,
            2,
            0,
            6,
        ),
        (["dir_r1_test/test_0.fq", "dir_r1_test/test_1.fq"], 1, 1, 0, 2),
        # 3 files, 2 samples
        (["test_0.fq", "other.fq", "test_1.fq"], 2, 2, 0, 2),
        # 2 paired files, no ONT parts
        (["test_R1.fq", "test_R2.fq"], 1, 2, 2, 0),
    ],
)
def test_paired_and_ont_files(
    runner,
    generate_fastq,
    mock_file_upload,
    mock_sample_get,
    mocked_creds_path,
    upload_mocks,
    files,
    n_samples_uploaded,
    n_files_uploaded,
    n_paired_files,
    n_ont_files,
):
    files = [generate_fastq(x) for x in files]

    args = ["--api-key", "01234567890123456789012345678901", "upload"] + files
    result = runner.invoke(Cli, args, catch_exceptions=False, input="\n\n")
    assert mock_file_upload.call_count == n_files_uploaded
    assert mock_sample_get.call_count == n_samples_uploaded
    assert result.exit_code == 0

    paired_files_prompt = "It appears there are {} paired files (of {} total)".format(
        n_paired_files, len(files)
    )
    if n_paired_files > 0:
        assert paired_files_prompt in result.output
    else:
        assert paired_files_prompt not in result.output

    ont_prompt = "Would you like to merge files by sample?"
    if n_ont_files > 0:
        assert ont_prompt in result.output
    else:
        assert ont_prompt not in result.output


def test_download_samples_without_prompt(runner, api_data, mocked_creds_file):
    with runner.isolated_filesystem():
        result = runner.invoke(Cli, ["download", "samples", "output", "--no-prompt"])
        assert result.exit_code == 0
        assert len(os.listdir("output")) > 0


def test_download_samples_with_prompt_confirmation(runner, api_data, mocked_creds_file):
    with runner.isolated_filesystem():
        result = runner.invoke(Cli, ["download", "samples", "output", "--prompt"], input="y\n")
        assert result.exit_code == 0
        assert len(os.listdir("output")) > 0


def test_download_samples_with_prompt_abort(runner, api_data, mocked_creds_file):
    with runner.isolated_filesystem():
        result = runner.invoke(Cli, ["download", "samples", "output", "--prompt"], input="n\n")
        assert result.exit_code == 0
        assert not os.path.exists("output")
