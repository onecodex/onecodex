# -*- coding: utf-8 -*-
from collections import OrderedDict
from io import BytesIO
from mock import patch
import pytest
from requests_toolbelt import MultipartEncoder
from requests.exceptions import HTTPError

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


class FakeAPISession:
    # TODO: We should migrate this to use responses vs. a faked API session
    adapters = {}

    def post(self, url, **kwargs):
        resp = lambda: None  # noqa
        resp.json = lambda: {"code": 200}
        resp.status_code = 201 if "auth" in kwargs else 200
        resp.raise_for_status = lambda: None  # noqa

        if url.endswith("s3_confirm"):
            resp.status_code = 200
            resp.json = lambda: {"sample_id": "s3_confirm_sample_id"}  # noqa

        return resp

    def mount(self, url, adapter):
        self.adapters[url] = adapter


class FakeAssetsResource:
    @staticmethod
    def err_resp():
        resp = lambda: None  # noqa
        resp.status_code = 400
        resp.json = lambda: {}  # noqa
        raise HTTPError(response=resp)

    def init_multipart_upload(self):
        return {
            "callback_url": "/s3_confirm",
            "s3_bucket": "some_bucket",
            "file_id": "hey",
            "upload_aws_access_key_id": "key",
            "upload_aws_secret_access_key": "secret",
        }

    def confirm_upload(self, obj):
        if self.what_fails == "confirm":
            self.err_resp()

        assert "sample_id" in obj
        assert "upload_type" in obj

    class _client(object):
        _root_url = "http://localhost:3000"
        session = FakeAPISession()


class FakeSamplesResource:
    def __init__(self, what_fails=None):
        self.what_fails = what_fails

    @staticmethod
    def err_resp():
        resp = lambda: None  # noqa
        resp.status_code = 400
        resp.json = lambda: {}  # noqa
        raise HTTPError(response=resp)

    def init_multipart_upload(self, obj):
        if self.what_fails == "init_multipart":
            self.err_resp()

        return {
            "callback_url": "/s3_confirm",
            "s3_bucket": "some_bucket",
            "file_id": "hey",
            "paired_end_file_id": "ho",
            "upload_aws_access_key_id": "key",
            "upload_aws_secret_access_key": "secret",
        }

    def confirm_upload(self, obj):
        if self.what_fails == "confirm":
            self.err_resp()

        assert "sample_id" in obj
        assert "upload_type" in obj

    def cancel_upload(self, obj):
        if self.what_fails == "cancel":
            self.err_resp()

        assert "sample_id" in obj
        return {"success": True, "message": "Upload cancelled"}

    class _client(object):
        _root_url = "http://localhost:3000"
        session = FakeAPISession()


class FakeDocumentsResource:
    def __init__(self, what_fails=None):
        self.what_fails = what_fails

    @staticmethod
    def err_resp():
        resp = lambda: None  # noqa
        resp.status_code = 400
        resp.json = lambda: {}  # noqa
        raise HTTPError(response=resp)

    def init_multipart_upload(self):
        if self.what_fails == "init_multipart":
            self.err_resp()

        return {
            "callback_url": "/s3_confirm",
            "s3_bucket": "some_bucket",
            "file_id": "hey",
            "upload_aws_access_key_id": "key",
            "upload_aws_secret_access_key": "secret",
        }

    def confirm_upload(self, obj):
        if self.what_fails == "confirm":
            self.err_resp()

        assert "sample_id" in obj
        assert "upload_type" in obj

    class _client(object):
        _root_url = "http://localhost:3000"
        session = FakeAPISession()


@pytest.mark.parametrize(
    "files,n_uploads,fxi_calls,fxp_calls",
    [
        ("file.1000.fa", 1, 0, 1),
        ("file.5e10.fa", 1, 0, 1),
        (("file.3e10.R1.fa", "file.3e10.R2.fa"), 2, 1, 0),
        # (["file.1000.fa", "file.5e10.fa"], 2, 0, 2, 2),
        # ([("file.3e10.R1.fa", "file.3e10.R2.fa")], 1, 1, 0, 2),
        # ([("file.3e5.R1.fa", "file.3e5.R2.fa"), "file.1e5.fa"], 2, 1, 1, 3),
    ],
)
def test_upload_lots_of_files(files, n_uploads, fxi_calls, fxp_calls):
    with patch("boto3.session.Session"):
        fake_size = lambda filename: int(float(filename.split(".")[1]))  # noqa

        uso = "onecodex.lib.upload._upload_sequence_fileobj"
        udo = "onecodex.lib.upload._upload_document_fileobj"
        fxi = "onecodex.lib.files.PairedEndFiles"
        fxp = "onecodex.lib.files.FilePassthru"
        sz = "os.path.getsize"

        sample_resource = FakeSamplesResource()
        with patch(uso) as upload_sequence_fileobj, patch(fxi) as paired, patch(fxp) as passthru:
            with patch(sz, size_effect=fake_size):
                upload_sequence(files, sample_resource)

                assert upload_sequence_fileobj.call_count == n_uploads
                assert paired.call_count == fxi_calls
                assert passthru.call_count == fxp_calls

        # Need to patch from `upload.` instead of `files.` for that test
        with patch(udo) as upload_document_fileobj, patch(
            "onecodex.lib.upload.FilePassthru"
        ) as passthru:
            with patch(sz, size_effect=fake_size):
                files = files[0] if isinstance(files, tuple) else files
                upload_document(files, FakeDocumentsResource())

                assert (
                    upload_document_fileobj.call_count == n_uploads - 1
                    if isinstance(files, tuple)
                    else n_uploads
                )
                assert passthru.call_count == fxp_calls + fxi_calls


def test_upload_asset():
    with patch("boto3.session.Session"):
        file = "test_asset_file.fa"
        n_uploads = 1
        fake_size = 1000

        with patch("onecodex.lib.upload._upload_asset_fileobj") as upload_asset_fileobj, patch(
            "onecodex.lib.upload.FilePassthru"
        ) as passthru, patch("os.path.getsize", size_effect=fake_size):
            upload_asset(file, FakeAssetsResource())

            assert upload_asset_fileobj.call_count == n_uploads
            assert upload_asset_fileobj.call_args[-1] == {"name": None}
            assert passthru.call_count == 1


def test_upload_asset_with_name():
    with patch("boto3.session.Session"):
        file = "test_asset_file.fa"
        fake_size = 1000
        name = "user-friendly-name"

        with patch("onecodex.lib.upload._upload_asset_fileobj") as upload_asset_fileobj, patch(
            "onecodex.lib.upload.FilePassthru"
        ), patch("os.path.getsize", size_effect=fake_size):
            upload_asset(file, FakeAssetsResource(), name=name)
            assert upload_asset_fileobj.call_args[-1] == {"name": name}


def test_api_failures(caplog):
    with patch("boto3.session.Session"):
        files = ("tests/data/files/test_R1_L001.fq.gz", "tests/data/files/test_R2_L001.fq.gz")

        with pytest.raises(UploadException) as e:
            upload_sequence(files, FakeSamplesResource("init_multipart"))
        assert "Could not initialize upload" in str(e.value)


def test_unicode_filenames(caplog):
    with patch("boto3.session.Session"):
        file_list = [
            ("tests/data/files/François.fq", "Francois.fq"),
            ("tests/data/files/Málaga.fasta", "Malaga.fasta"),
            ("tests/data/files/Röö.fastq", "Roo.fastq"),
        ]

        # should raise if --coerce-ascii not passed
        sample_resource = FakeSamplesResource()
        for before, after in file_list:
            with pytest.raises(OneCodexException) as e:
                upload_sequence(before, sample_resource)
            assert "must be ascii" in str(e.value)

        # make sure log gets warnings when we rename files
        for before, _ in file_list:
            upload_sequence(before, sample_resource, coerce_ascii=True)

        for _, after in file_list:
            assert after in caplog.text


def test_single_end_files():
    with patch("boto3.session.Session") as b3:
        upload_sequence("tests/data/files/test_R1_L001.fq.gz", FakeSamplesResource())
        assert b3.call_count == 1


def test_paired_end_files():
    with patch("boto3.session.Session") as b3:
        upload_sequence(
            ("tests/data/files/test_R1_L001.fq.gz", "tests/data/files/test_R2_L001.fq.gz"),
            FakeSamplesResource(),
        )
        assert b3.call_count == 2


def test_upload_sequence_fileobj():
    with patch("boto3.session.Session") as b3:
        # upload succeeds via multipart
        file_obj = BytesIO(b">test\nACGT\n")
        _upload_sequence_fileobj(
            file_obj,
            "test.fa",
            FakeSamplesResource().init_multipart_upload({}),
            FakeSamplesResource(),
        )
        file_obj.close()
        assert b3.call_count == 1


def test_upload_document_fileobj():
    with patch("boto3.session.Session") as b3:
        file_obj = BytesIO(b"MY_SPOON_IS_TOO_BIG\n")
        _upload_document_fileobj(file_obj, "spoon.pdf", FakeDocumentsResource())
        file_obj.close()

        assert b3.call_count == 1

    with pytest.raises(UploadException) as e:
        file_obj = BytesIO(b"MY_SPOON_IS_TOO_BIG\n")
        _upload_document_fileobj(file_obj, "spoon.pdf", FakeDocumentsResource("init_multipart"))
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
    assert (
        _choose_boto3_chunksize(open("tests/data/files/test_R1_L001.fq.gz", "r")) == 25 * 1024**2
    )
