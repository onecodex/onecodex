import os
import os.path

import pytest

from onecodex.exceptions import OneCodexException
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


def test_download_samples_does_not_overwrite_existing_files(runner, ocx, api_data):
    with runner.isolated_filesystem():
        download_samples(ocx, "output")

        with pytest.raises(OneCodexException):
            download_samples(ocx, "output")


def test_download_samples_outdir_exists(runner, ocx, api_data):
    with runner.isolated_filesystem():
        os.mkdir("output")
        download_samples(ocx, "output")
