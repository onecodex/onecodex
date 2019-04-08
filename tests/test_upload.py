# -*- coding: utf-8 -*-
from collections import OrderedDict
from hashlib import sha256
from io import BytesIO
import logging
from mock import patch
import os
import pytest
from requests_toolbelt import MultipartEncoder
from requests.exceptions import HTTPError

from onecodex.exceptions import OneCodexException, UploadException
from onecodex.lib.upload import (
    _file_stats,
    connectivity_error,
    FASTXInterleave,
    FilePassthru,
    interleaved_filename,
    upload_document,
    upload_document_fileobj,
    upload_sequence,
    upload_sequence_fileobj,
)


@pytest.mark.parametrize(
    "files,filename",
    [
        (("test_R1.fastq", "test_R2.fastq"), "test.fastq"),
        (("test_R1_001.fastq", "test_R2_001.fastq"), "test_001.fastq"),
    ],
)
@pytest.mark.filterwarnings("error")
def test_interleaved_filenames(files, filename):
    assert interleaved_filename(files) == filename

    with pytest.raises(OneCodexException) as e:
        interleaved_filename(files[0])
    assert "without a tuple" in str(e.value)

    with pytest.raises(UserWarning) as e:
        interleaved_filename(tuple([x.replace("_R", "") for x in files]))
    assert "do not match" in str(e.value)


def test_file_stats():
    pair = ("tests/data/files/test_R1_L001.fq.gz", "tests/data/files/test_R2_L001.fq.gz")
    stats = _file_stats(pair)
    assert stats == ("test_L001.fq", 2420, "fastq")

    with pytest.raises(UploadException) as e:
        _file_stats("tests/data/api/complete_analyses.json")
    assert "extension must be one of" in str(e.value)

    stats = _file_stats("tests/data/api/complete_analyses.json", enforce_fastx=False)
    assert stats == ("complete_analyses.json", 683, None)


def test_fastxinterleave():
    pair = ("tests/data/files/test_R1_L001.fq.gz", "tests/data/files/test_R2_L001.fq.gz")
    wrapper = FASTXInterleave(pair, 2420, "fastq")
    assert (
        sha256(wrapper.read()).hexdigest()
        == "47a98ffeb62e53af916a7afd75e13d1d189d683d0e796c2eef9abb501869de0f"
    )

    pair = ("tests/data/files/test_R1_L001.fq.bz2", "tests/data/files/test_R2_L001.fq.bz2")
    wrapper = FASTXInterleave(pair, 2420, "fastq")
    assert (
        sha256(wrapper.read()).hexdigest()
        == "47a98ffeb62e53af916a7afd75e13d1d189d683d0e796c2eef9abb501869de0f"
    )

    with pytest.raises(OneCodexException) as e:
        wrapper = FASTXInterleave(pair, 2420, "fasta")
    assert "currently unsupported" in str(e.value)

    with pytest.raises(OneCodexException) as e:
        wrapper = FASTXInterleave(pair, 2420, "stockholm")
    assert "must be one of" in str(e.value)


class FakeSamplesResource:
    def __init__(self, what_fails=None):
        self.what_fails = what_fails

    @staticmethod
    def err_resp():
        resp = lambda: None  # noqa
        resp.status_code = 400
        resp.json = lambda: {}  # noqa
        raise HTTPError(response=resp)

    def init_upload(self, obj):
        if self.what_fails == "init":
            self.err_resp()

        assert "filename" in obj
        assert "size" in obj
        assert "upload_type" in obj
        return {
            "upload_url": "http://localhost:3005/fastx_proxy",
            "sample_id": "sample_uuid_here",
            "additional_fields": {"status_url": "http://localhost:3005/fastx_proxy/errors"},
        }

    def init_multipart_upload(self, obj):
        if self.what_fails == "init_multipart":
            self.err_resp()

        return {
            "callback_url": "/s3_confirm",
            "s3_bucket": "",
            "file_id": "",
            "upload_aws_access_key_id": "",
            "upload_aws_secret_access_key": "",
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
            "s3_bucket": "",
            "file_id": "",
            "upload_aws_access_key_id": "",
            "upload_aws_secret_access_key": "",
        }

    def confirm_upload(self, obj):
        if self.what_fails == "confirm":
            self.err_resp()

        assert "sample_id" in obj
        assert "upload_type" in obj

    class _client(object):
        _root_url = "http://localhost:3000"


class FakeSession:
    def post(self, url, **kwargs):
        if "data" in kwargs:
            kwargs["data"].read()  # So multipart uploader will close properly
        resp = lambda: None  # noqa
        resp.json = lambda: {}  # noqa
        resp.status_code = 201 if "auth" in kwargs else 200
        return resp


class FakeSessionProxyFails:
    def __init__(self, failure_code, no_msg=False):
        self.failure_code = failure_code
        self.no_msg = no_msg

    # proxy upload will fail and fall through to S3
    def post(self, url, **kwargs):
        if "data" in kwargs:
            kwargs["data"].read()  # So multipart uploader will close properly
        resp = lambda: None  # noqa

        if url.endswith("fastx_proxy"):
            resp.status_code = self.failure_code
            if self.no_msg:
                resp.json = lambda: {}  # noqa
            else:
                resp.json = lambda: {"message": "FASTX Proxy Error Message"}  # noqa
        elif url.endswith("s3_confirm"):
            resp.status_code = 200
            resp.json = lambda: {"sample_id": "s3_confirm_sample_id"}  # noqa
        elif "auth" in kwargs:
            resp.status_code = 201
        else:
            resp.status_code = 200

        # support for callback for more info on fastx-proxy errors
        if url.endswith("errors"):
            if self.no_msg:
                resp.status_code = 500
            else:
                resp.json = lambda: {"message": "Additional error message"}  # noqa

        return resp


@pytest.mark.parametrize(
    "file_list,nfiles,fxi_calls,fxp_calls,size_calls",
    [
        (["file.1000.fa", "file.5e10.fa"], 2, 0, 2, 2),
        ([("file.3e10.R1.fa", "file.3e10.R2.fa")], 1, 1, 0, 2),
        ([("file.3e5.R1.fa", "file.3e5.R2.fa"), "file.1e5.fa"], 2, 1, 1, 3),
    ],
)
def test_upload_lots_of_files(file_list, nfiles, fxi_calls, fxp_calls, size_calls):
    fake_size = lambda filename: int(float(filename.split(".")[1]))  # noqa

    uso = "onecodex.lib.upload.upload_sequence_fileobj"
    udo = "onecodex.lib.upload.upload_document_fileobj"
    fxi = "onecodex.lib.upload.FASTXInterleave"
    fxp = "onecodex.lib.upload.FilePassthru"
    sz = "onecodex.lib.upload.os.path.getsize"

    with patch(uso) as upload_sequence_fileobj, patch(fxi) as interleave, patch(fxp) as passthru:
        with patch(sz, size_effect=fake_size) as size:
            upload_sequence(file_list, FakeSession(), FakeSamplesResource(), threads=1)

            assert upload_sequence_fileobj.call_count == nfiles
            assert interleave.call_count == fxi_calls
            assert passthru.call_count == fxp_calls
            assert size.call_count == size_calls

    with patch(udo) as upload_document_fileobj, patch(fxp) as passthru:
        with patch(sz, size_effect=fake_size) as size:
            file_list = [(x[0] if isinstance(x, tuple) else x) for x in file_list]
            upload_document(file_list, FakeSession(), FakeDocumentsResource(), threads=1)

            assert upload_document_fileobj.call_count == nfiles
            assert passthru.call_count == fxp_calls + fxi_calls
            assert size.call_count == fxp_calls + fxi_calls


@pytest.mark.parametrize(
    "file_list,nfiles,fxi_calls,fxp_calls,size_calls",
    [
        (["file.1000.fa", "file.5e10.fa"], 2, 0, 2, 2),
        ([("file.3e10.R1.fa", "file.3e10.R2.fa")], 1, 1, 0, 2),
        ([("file.3e5.R1.fa", "file.3e5.R2.fa"), "file.1e5.fa"], 2, 1, 1, 3),
    ],
)
def test_upload_multithread(file_list, nfiles, fxi_calls, fxp_calls, size_calls):
    fake_size = lambda filename: int(float(filename.split(".")[1]))  # noqa

    uso = "onecodex.lib.upload.upload_sequence_fileobj"
    udo = "onecodex.lib.upload.upload_document_fileobj"
    fxi = "onecodex.lib.upload.FASTXInterleave"
    fxp = "onecodex.lib.upload.FilePassthru"
    sz = "onecodex.lib.upload.os.path.getsize"

    with patch(uso) as upload_sequence_fileobj, patch(fxi) as interleave, patch(fxp) as passthru:
        with patch(sz, size_effect=fake_size) as size:
            upload_sequence(file_list, FakeSession(), FakeSamplesResource(), threads=4)

            assert upload_sequence_fileobj.call_count == nfiles
            assert interleave.call_count == fxi_calls
            assert passthru.call_count == fxp_calls
            assert size.call_count == size_calls

    with patch(udo) as upload_document_fileobj, patch(fxp) as passthru:
        with patch(sz, size_effect=fake_size) as size:
            file_list = [(x[0] if isinstance(x, tuple) else x) for x in file_list]
            upload_document(file_list, FakeSession(), FakeDocumentsResource(), threads=4)

            assert upload_document_fileobj.call_count == nfiles
            assert passthru.call_count == fxp_calls + fxi_calls
            assert size.call_count == fxp_calls + fxi_calls


def test_api_failures(caplog):
    file_list = [
        ("tests/data/files/test_R1_L001.fq.gz", "tests/data/files/test_R2_L001.fq.gz"),
        ("tests/data/files/test_R1_L001.fq.bz2", "tests/data/files/test_R2_L001.fq.bz2"),
        "tests/data/files/test.fa",
    ]

    with pytest.raises(UploadException) as e:
        upload_sequence(file_list, FakeSession(), FakeSamplesResource("init"), threads=1)
    assert "Could not initialize upload" in str(e.value)

    upload_sequence(
        file_list, FakeSession(), FakeSamplesResource("confirm"), threads=1, log=logging
    )
    assert "Callback could not be completed" in caplog.text

    upload_sequence(
        file_list,
        FakeSessionProxyFails(400, no_msg=True),
        FakeSamplesResource(),
        threads=1,
        log=logging,
    )
    assert "File could not be uploaded" in caplog.text

    with pytest.raises(UploadException) as e:
        connectivity_error("filename")
    assert "experiencing connectivity" in str(e.value)


def test_unicode_filenames(caplog):
    file_list = [
        (u"tests/data/files/François.fq", "Francois.fq"),
        (u"tests/data/files/Málaga.fasta", "Malaga.fasta"),
        (u"tests/data/files/Röö.fastq", "Roo.fastq"),
    ]

    # should raise if --coerce-ascii not passed
    for before, after in file_list:
        with pytest.raises(OneCodexException) as e:
            upload_sequence([before], FakeSession(), FakeSamplesResource())
        assert "must be ascii" in str(e.value)

    # make sure log gets warnings when we rename files
    upload_sequence(
        [x[0] for x in file_list],
        FakeSession(),
        FakeSamplesResource(),
        log=logging,
        coerce_ascii=True,
    )

    for _, after in file_list:
        assert after in caplog.text


def test_single_end_files():
    upload_sequence(
        ["tests/data/files/test_R1_L001.fq.gz", "tests/data/files/test_R2_L001.fq.gz"],
        FakeSession(),
        FakeSamplesResource(),
        threads=1,
    )


def test_paired_end_files():
    upload_sequence(
        [("tests/data/files/test_R1_L001.fq.gz", "tests/data/files/test_R2_L001.fq.gz")],
        FakeSession(),
        FakeSamplesResource(),
        threads=1,
    )


def test_upload_sequence_fileobj():
    # upload succeeds via fastx-proxy
    file_obj = BytesIO(b">test\nACGT\n")
    init_upload_fields = FakeSamplesResource().init_upload(["filename", "size", "upload_type"])
    upload_sequence_fileobj(
        file_obj, "test.fa", init_upload_fields, {}, FakeSession(), FakeSamplesResource()
    )
    file_obj.close()

    # upload via fastx-proxy fails via 500 and must fall through to S3
    with patch("boto3.client") as b3:
        file_obj = BytesIO(b">test\nACGT\n")
        init_upload_fields = FakeSamplesResource().init_upload(["filename", "size", "upload_type"])
        ret_sample_id = upload_sequence_fileobj(
            file_obj,
            "test.fa",
            init_upload_fields,
            {},
            FakeSessionProxyFails(500),
            FakeSamplesResource(),
        )
        file_obj.close()

        assert b3.call_count == 1
        assert ret_sample_id == "s3_confirm_sample_id"

    # upload via fastx-proxy fails via 400 (a legitimate, probably validation error), and we're able
    # to contact the API to get a more detailed error message
    with pytest.raises(UploadException) as e:
        file_obj = BytesIO(b">test\nACGT\n")
        init_upload_fields = FakeSamplesResource().init_upload(["filename", "size", "upload_type"])
        upload_sequence_fileobj(
            file_obj,
            "test.fa",
            init_upload_fields,
            {},
            FakeSessionProxyFails(400),
            FakeSamplesResource(),
        )
        file_obj.close()
    assert "Additional error message" in str(e.value)


def test_upload_document_fileobj():
    with patch("boto3.client") as b3:
        file_obj = BytesIO(b"MY_SPOON_IS_TOO_BIG\n")
        upload_document_fileobj(file_obj, "spoon.pdf", FakeSession(), FakeDocumentsResource())
        file_obj.close()

        assert b3.call_count == 1

    with pytest.raises(UploadException) as e:
        file_obj = BytesIO(b"MY_SPOON_IS_TOO_BIG\n")
        upload_document_fileobj(
            file_obj, "spoon.pdf", FakeSession(), FakeDocumentsResource("init_multipart")
        )
        file_obj.close()
    assert "Could not initialize" in str(e.value)


def test_multipart_encoder_passthru():
    wrapper = FilePassthru(
        "tests/data/files/test_R1_L001.fq.bz2",
        os.path.getsize("tests/data/files/test_R1_L001.fq.bz2"),
    )
    wrapper_len = len(wrapper.read())
    wrapper.seek(0)

    assert wrapper_len == wrapper._fsize

    multipart_fields = OrderedDict()
    multipart_fields["file"] = ("fakefile", wrapper, "application/x-gzip")
    encoder = MultipartEncoder(multipart_fields)
    MAGIC_HEADER_LEN = 178
    wrapper.seek(0)
    assert len(encoder.read()) - MAGIC_HEADER_LEN == wrapper_len


def test_multipart_encoder_interleave():
    pair = ("tests/data/files/test_R1_L001.fq.gz", "tests/data/files/test_R2_L001.fq.gz")
    fname, fsize, fformat = _file_stats(pair)
    wrapper = FASTXInterleave(pair, fsize, fformat)
    wrappertext = wrapper.read()
    wrapper.seek(0)

    assert len(wrappertext) == fsize

    multipart_fields = OrderedDict()
    multipart_fields["file"] = ("fakefile", wrapper, "text/plain")
    encoder = MultipartEncoder(multipart_fields)
    MAGIC_HEADER_LEN = 170  # shorter because of text/plain mime-type
    encodertext = encoder.read()
    assert len(encodertext) - MAGIC_HEADER_LEN == len(wrappertext)
