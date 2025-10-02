# -*- coding: utf-8 -*-
from collections import OrderedDict
from io import BytesIO
from mock import patch
import pytest
import responses
from requests_toolbelt import MultipartEncoder

from onecodex.exceptions import OneCodexException, UploadException
from onecodex.lib.upload import (
    _choose_boto3_chunksize,
    FilePassthru,
    upload_document,
    _upload_document_fileobj,
    upload_asset,
    upload_sequence,
    _upload_sequence_fileobj,
)


@pytest.fixture
def mock_api_responses():
    """Setup responses for successful API calls."""
    with responses.RequestsMock(assert_all_requests_are_fired=False) as rsps:
        # Samples endpoints
        rsps.add(
            responses.POST,
            "http://localhost:3000/api/v1/samples/init_multipart_upload",
            json={
                "callback_url": "/s3_confirm",
                "s3_bucket": "some_bucket",
                "file_id": "hey",
                "paired_end_file_id": "ho",
                "upload_aws_access_key_id": "key",
                "upload_aws_secret_access_key": "secret",
                "sample_id": "test_sample_id",
            },
            status=200,
        )

        rsps.add(
            responses.POST,
            "http://localhost:3000/s3_confirm",
            json={"sample_id": "abcdef1234567890"},
            status=200,
        )

        rsps.add(
            responses.GET,
            "http://localhost:3000/api/v1/samples/abcdef1234567890",
            json={
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
                "status": "available",
            },
            status=200,
        )

        rsps.add(
            responses.POST,
            "http://localhost:3000/api/v1/samples/cancel_upload",
            json={"success": True, "message": "Upload cancelled"},
            status=200,
        )

        # Documents endpoints
        rsps.add(
            responses.POST,
            "http://localhost:3000/api/v1/documents/init_multipart_upload",
            json={
                "callback_url": "/s3_confirm",
                "s3_bucket": "some_bucket",
                "file_id": "hey",
                "upload_aws_access_key_id": "key",
                "upload_aws_secret_access_key": "secret",
            },
            status=200,
        )

        # Assets endpoints
        rsps.add(
            responses.POST,
            "http://localhost:3000/api/v1/assets/init_multipart_upload",
            json={
                "callback_url": "/s3_confirm",
                "s3_bucket": "some_bucket",
                "file_id": "hey",
                "upload_aws_access_key_id": "key",
                "upload_aws_secret_access_key": "secret",
            },
            status=200,
        )

        yield rsps


@pytest.fixture
def mock_api_failures():
    """Setup responses for API failure scenarios."""
    with responses.RequestsMock(assert_all_requests_are_fired=False) as rsps:
        rsps.add(
            responses.POST,
            "http://localhost:3000/api/v1/samples/init_multipart_upload",
            json={},
            status=500,
        )
        rsps.add(
            responses.POST,
            "http://localhost:3000/api/v1/documents/init_multipart_upload",
            json={},
            status=500,
        )

        yield rsps


@pytest.mark.parametrize(
    "files,n_uploads,fxi_calls,fxp_calls",
    [
        ("file.1000.fa", 1, 0, 1),
        ("file.5e10.fa", 1, 0, 1),
        (("file.3e10.R1.fa", "file.3e10.R2.fa"), 2, 1, 0),
    ],
)
def test_upload_lots_of_files(files, n_uploads, fxi_calls, fxp_calls, ocx, mock_api_responses):
    with patch("boto3.session.Session"):
        fake_size = lambda filename: int(float(filename.split(".")[1]))  # noqa

        uso = "onecodex.lib.upload._upload_sequence_fileobj"
        udo = "onecodex.lib.upload._upload_document_fileobj"
        fxi = "onecodex.lib.files.PairedEndFiles"
        fxp = "onecodex.lib.files.FilePassthru"
        sz = "os.path.getsize"

        with patch(uso) as upload_sequence_fileobj, patch(fxi) as paired, patch(fxp) as passthru:
            # Configure the mocks to return proper filenames for JSON serialization
            if isinstance(files, tuple):
                # Paired-end files
                paired_mock = paired.return_value
                paired_mock.r1.filename = files[0]
                paired_mock.r2.filename = files[1]
            else:
                # Single file
                passthru_mock = passthru.return_value
                passthru_mock.filename = files

            with patch(sz, side_effect=fake_size):
                upload_sequence(files, ocx.Samples)

                assert upload_sequence_fileobj.call_count == n_uploads
                assert paired.call_count == fxi_calls
                assert passthru.call_count == fxp_calls

        # Need to patch from `upload.` instead of `files.` for that test
        with (
            patch(udo) as upload_document_fileobj,
            patch("onecodex.lib.upload.FilePassthru") as passthru,
        ):
            # Configure filename for document upload
            passthru_mock = passthru.return_value
            passthru_mock.filename = files[0] if isinstance(files, tuple) else files

            with patch(sz, side_effect=fake_size):
                files = files[0] if isinstance(files, tuple) else files
                upload_document(files, ocx.Documents)

                assert (
                    upload_document_fileobj.call_count == n_uploads - 1
                    if isinstance(files, tuple)
                    else n_uploads
                )
                assert passthru.call_count == fxp_calls + fxi_calls


def test_upload_asset(ocx, mock_api_responses):
    with patch("boto3.session.Session"):
        file = "test_asset_file.fa"
        n_uploads = 1
        fake_size = 1000

        with (
            patch("onecodex.lib.upload._upload_asset_fileobj") as upload_asset_fileobj,
            patch("onecodex.lib.upload.FilePassthru") as passthru,
            patch("os.path.getsize", side_effect=lambda x: fake_size),
        ):
            # Configure the mock to have a proper filename
            passthru_mock = passthru.return_value
            passthru_mock.filename = file

            upload_asset(file, ocx.Assets)

            assert upload_asset_fileobj.call_count == n_uploads
            assert upload_asset_fileobj.call_args[-1] == {"name": None}
            assert passthru.call_count == 1


def test_upload_asset_with_name(ocx, mock_api_responses):
    with patch("boto3.session.Session"):
        file = "test_asset_file.fa"
        fake_size = 1000
        name = "user-friendly-name"

        with (
            patch("onecodex.lib.upload._upload_asset_fileobj") as upload_asset_fileobj,
            patch("onecodex.lib.upload.FilePassthru") as passthru,
            patch("os.path.getsize", side_effect=lambda x: fake_size),
        ):
            # Configure the mock to have a proper filename
            passthru_mock = passthru.return_value
            passthru_mock.filename = file

            upload_asset(file, ocx.Assets, name=name)
            assert upload_asset_fileobj.call_args[-1] == {"name": name}


def test_api_failures(caplog, ocx, mock_api_failures):
    with (
        patch("onecodex.lib.upload.get_file_wrapper") as mock_wrapper,
        patch("boto3.session.Session"),
    ):
        files = ("tests/data/files/test_R1_L001.fq.gz", "tests/data/files/test_R2_L001.fq.gz")

        # Configure mock paired-end file object with required properties
        mock_file_obj = mock_wrapper.return_value
        mock_file_obj.r1.filename = files[0]
        mock_file_obj.r2.filename = files[1]
        mock_file_obj._fsize = 1000  # Add size property to avoid comparison errors

        with pytest.raises(UploadException) as e:
            upload_sequence(files, ocx.Samples)
        assert "Could not initialize upload" in str(e.value)


def test_unicode_filenames(caplog, ocx, mock_api_responses):
    # Don't mock get_file_wrapper so that ASCII validation actually happens
    with patch("boto3.session.Session"):
        file_list = [
            ("tests/data/files/François.fq", "Francois.fq"),
            ("tests/data/files/Málaga.fasta", "Malaga.fasta"),
            ("tests/data/files/Röö.fastq", "Roo.fastq"),
        ]

        # should raise if --coerce-ascii not passed
        for before, after in file_list:
            with pytest.raises(OneCodexException) as e:
                upload_sequence(before, ocx.Samples)
            assert "must be ascii" in str(e.value)

        # make sure log gets warnings when we rename files
        for before, _ in file_list:
            upload_sequence(before, ocx.Samples, coerce_ascii=True)

        for _, after in file_list:
            assert after in caplog.text


def test_single_end_files(ocx, mock_api_responses):
    with (
        patch("onecodex.lib.upload.get_file_wrapper") as mock_wrapper,
        patch("boto3.session.Session") as b3,
    ):
        # Configure mock file object
        mock_file_obj = mock_wrapper.return_value
        mock_file_obj.filename = "test_R1_L001.fq.gz"
        mock_file_obj._fsize = 1000  # Add size property to avoid comparison errors

        upload_sequence("tests/data/files/test_R1_L001.fq.gz", ocx.Samples)
        assert b3.call_count == 1


def test_single_end_files_upload_method(ocx, mock_api_responses):
    with (
        patch("onecodex.lib.upload.get_file_wrapper") as mock_wrapper,
        patch("boto3.session.Session") as b3,
    ):
        # Configure mock file object
        mock_file_obj = mock_wrapper.return_value
        mock_file_obj.filename = "test_R1_L001.fq.gz"
        mock_file_obj._fsize = 1000  # Add size property to avoid comparison errors

        ocx.Samples.upload("tests/data/files/test_R1_L001.fq.gz")
        assert b3.call_count == 1


def test_paired_end_files(ocx, mock_api_responses):
    with (
        patch("onecodex.lib.upload.get_file_wrapper") as mock_wrapper,
        patch("boto3.session.Session") as b3,
    ):
        # Configure mock paired-end file object
        mock_file_obj = mock_wrapper.return_value
        mock_file_obj.r1.filename = "test_R1_L001.fq.gz"
        mock_file_obj.r2.filename = "test_R2_L001.fq.gz"
        # Set _fsize on the individual r1 and r2 objects for _choose_boto3_chunksize
        mock_file_obj.r1._fsize = 1000
        mock_file_obj.r2._fsize = 1000
        mock_file_obj._fsize = 1000  # Also set on parent for safety

        upload_sequence(
            ("tests/data/files/test_R1_L001.fq.gz", "tests/data/files/test_R2_L001.fq.gz"),
            ocx.Samples,
        )
        assert b3.call_count == 2


def test_upload_sequence_fileobj(ocx, mock_api_responses):
    with patch("boto3.session.Session") as b3:
        # upload succeeds via multipart
        file_obj = BytesIO(b">test\nACGT\n")
        fields = {
            "callback_url": "/s3_confirm",
            "s3_bucket": "some_bucket",
            "file_id": "hey",
            "upload_aws_access_key_id": "key",
            "upload_aws_secret_access_key": "secret",
        }
        _upload_sequence_fileobj(
            file_obj,
            "test.fa",
            fields,
            ocx.Samples,
        )
        file_obj.close()
        assert b3.call_count == 1


def test_upload_document_fileobj(ocx, mock_api_responses):
    # Test successful upload first (with success responses active)
    with patch("boto3.session.Session") as b3:
        file_obj = BytesIO(b"MY_SPOON_IS_TOO_BIG\n")
        _upload_document_fileobj(file_obj, "spoon.pdf", ocx.Documents)
        file_obj.close()
        assert b3.call_count == 1


def test_upload_document_fileobj_failure(ocx, mock_api_failures):
    # Test failure case separately with failure responses
    with patch("boto3.session.Session"):
        with pytest.raises(UploadException) as e:
            file_obj = BytesIO(b"MY_SPOON_IS_TOO_BIG\n")
            _upload_document_fileobj(file_obj, "spoon.pdf", ocx.Documents)
            file_obj.close()
        assert "Could not initialize" in str(e.value)


def test_multipart_encoder_passthru():
    wrapper = FilePassthru("tests/data/files/test_R1_L001.fq.bz2")
    wrapper_len = len(wrapper.read())
    wrapper.seek(0)

    assert wrapper_len == wrapper._fsize

    multipart_fields = OrderedDict()
    multipart_fields["file"] = ("fakefile", wrapper, "application/x-gzip")
    encoder = MultipartEncoder(multipart_fields)
    MAGIC_HEADER_LEN = 178
    wrapper.seek(0)
    assert len(encoder.read()) - MAGIC_HEADER_LEN == wrapper_len


def test_boto3_chunksize():
    # test file that is too large to upload
    wrapper = FilePassthru("tests/data/files/test_R1_L001.fq.gz")
    wrapper._fsize = 1024**4  # 1 TB

    with pytest.raises(OneCodexException) as e:
        _choose_boto3_chunksize(wrapper)
    assert "too large to upload" in str(e.value)

    # test file with no known size
    assert _choose_boto3_chunksize(open("tests/data/files/test_R1_L001.fq.gz", "r")) == 25 * 1024**2
