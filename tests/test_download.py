import os
import os.path

import pytest

from onecodex.exceptions import OneCodexException
from onecodex.lib.download import (
    download_samples,
    filter_samples_by_tags,
    get_project,
    get_download_filename,
)


class MockApi(object):
    def __init__(self, project_ids, project_names, tag_names):
        self.Projects = MockProjects(project_ids, project_names)
        self.Tags = MockTags(tag_names)


class MockProjects(object):
    def __init__(self, project_ids, project_names):
        self._project_ids = project_ids
        self._project_names = project_names

    def get(self, project_id):
        if project_id in self._project_ids:
            return project_id
        return None

    def where(self, name=None):
        if name and name in self._project_names:
            return [name]
        return []


class MockTags(object):
    def __init__(self, tag_names):
        self._tag_names = tag_names

    def where(self, name=None):
        if name and name in self._tag_names:
            return [name]
        return []


class MockSamples(object):
    def __init__(self, id="", filename="", tag_names=None):
        self.id = id
        self.filename = filename
        self.tags = tag_names if tag_names else []


@pytest.fixture
def mock_api():
    return MockApi(
        project_ids=["id1", "id2"],
        project_names=["proj1", "proj2"],
        tag_names=["tag1", "tag2", "tag3"],
    )


def test_download_samples(runner, ocx, api_data):
    with runner.isolated_filesystem():
        filepaths = download_samples(ocx, "output")
        for filepath in filepaths:
            assert os.path.exists(filepath)


def test_download_samples_with_progressbar(runner, ocx, api_data):
    with runner.isolated_filesystem():
        filepaths = download_samples(ocx, "output", progressbar=True)
        for filepath in filepaths:
            assert os.path.exists(filepath)


@pytest.mark.filterwarnings("ignore:Skipping download of sample.*")
def test_download_samples_does_not_overwrite_existing_files(runner, ocx, api_data):
    with runner.isolated_filesystem():
        filepaths1 = download_samples(ocx, "output")
        assert len(filepaths1) > 0

        filepaths2 = download_samples(ocx, "output")
        assert len(filepaths2) == 0


def test_download_samples_outdir_exists(runner, ocx, api_data):
    with runner.isolated_filesystem():
        os.mkdir("output")
        download_samples(ocx, "output")


def test_get_project_by_id(mock_api):
    project = get_project(mock_api, "id2")
    assert project == "id2"


def test_get_project_by_name(mock_api):
    project = get_project(mock_api, "proj1")
    assert project == "proj1"


def test_get_project_does_not_exist(mock_api):
    with pytest.raises(OneCodexException):
        get_project(mock_api, "id3")


@pytest.mark.filterwarnings("ignore:No tag found.*")
def test_filter_samples_by_tags(mock_api):
    samples = [
        MockSamples(),
        MockSamples(tag_names=["tag1"]),
        MockSamples(tag_names=["tag2", "tag3"]),
        MockSamples(tag_names=["tag3"]),
    ]

    filtered_samples = filter_samples_by_tags(mock_api, samples, ["tag1", "tag2", "tag4"])

    assert filtered_samples == [samples[1], samples[2]]

    # No tags found
    filtered_samples = filter_samples_by_tags(mock_api, samples, ["tag4", "tag5"])

    assert filtered_samples == samples


@pytest.mark.parametrize(
    "sample_id,sample_filename,output_filename",
    [
        ("abc123", "foo", "foo_abc123"),
        ("abc123", "foo.fastq", "foo_abc123.fastq"),
        ("abc123", "foo.fastq.gz", "foo_abc123.fastq.gz"),
    ],
)
def test_get_download_filename(sample_id, sample_filename, output_filename):
    sample = MockSamples(id=sample_id, filename=sample_filename)
    filename = get_download_filename(sample)

    assert filename == output_filename
