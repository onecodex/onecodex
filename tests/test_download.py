import os
import os.path

import pytest

from onecodex.lib.download import download_samples


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
