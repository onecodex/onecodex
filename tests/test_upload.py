# -*- coding: utf-8 -*-
from collections import OrderedDict
from io import BytesIO
import logging
from mock import patch
import os
import pytest
from requests_toolbelt import MultipartEncoder

from onecodex.exceptions import OneCodexException
from onecodex.lib.upload import (upload, upload_fileobj, interleaved_filename, FASTXPassthru,
                                 FASTXInterleave, _file_stats)


@pytest.mark.parametrize('files,filename', [
    (('test_R1.fastq', 'test_R2.fastq'), 'test.fastq'),
    (('test_R1_001.fastq', 'test_R2_001.fastq'), 'test_001.fastq')
])
def test_interleaved_filenames(files, filename):
    assert interleaved_filename(files) == filename


@pytest.mark.parametrize('file_list,nfiles,fxi_calls,fxp_calls,size_calls', [
    (['file.1000.fa', 'file.5e10.fa'], 2, 0, 2, 2),
    ([('file.3e10.R1.fa', 'file.3e10.R2.fa')], 1, 1, 0, 2),
    ([('file.3e5.R1.fa', 'file.3e5.R2.fa'), 'file.1e5.fa'], 2, 1, 1, 3),
])
def test_upload_lots_of_files(file_list, nfiles, fxi_calls, fxp_calls, size_calls):
    session, samples_resource = None, None

    fake_size = lambda filename: int(float(filename.split('.')[1]))  # noqa

    ufo = 'onecodex.lib.upload.upload_fileobj'
    fxi = 'onecodex.lib.upload.FASTXInterleave'
    fxp = 'onecodex.lib.upload.FASTXPassthru'
    sz = 'onecodex.lib.upload.os.path.getsize'

    with patch(ufo) as upload_fileobj, patch(fxi) as interleave, patch(fxp) as passthru:
        with patch(sz, size_effect=fake_size) as size:
            upload(file_list, session, samples_resource, threads=1)

            assert upload_fileobj.call_count == nfiles
            assert interleave.call_count == fxi_calls
            assert passthru.call_count == fxp_calls
            assert size.call_count == size_calls


class FakeSamplesResource():
    def init_upload(self, obj):
        assert 'filename' in obj
        assert 'size' in obj
        assert 'upload_type' in obj
        return {
            'upload_url': '',
            'sample_id': '',
            'additional_fields': {},
        }

    def init_multipart_upload(self):
        return {
            'callback_url': '',
            's3_bucket': '',
            'file_id': '',
            'upload_aws_access_key_id': '',
            'upload_aws_secret_access_key': '',
        }

    def confirm_upload(self, obj):
        assert 'sample_id' in obj
        assert 'upload_type' in obj


class FakeSession():
    def post(self, url, **kwargs):
        if 'data' in kwargs:
            kwargs['data'].read()  # So multipart uploader will close properly
        resp = lambda: None  # noqa
        resp.status_code = 201 if 'auth' in kwargs else 200
        return resp


def test_unicode_filenames(caplog):
    file_list = [
        (u'tests/data/files/François.fq', 'Francois.fq'),
        (u'tests/data/files/Málaga.fasta', 'Malaga.fasta'),
        (u'tests/data/files/Röö.fastq', 'Roo.fastq')
    ]

    # should raise if --coerce-ascii not passed
    for before, after in file_list:
        with pytest.raises(OneCodexException) as e:
            upload([before], FakeSession(), FakeSamplesResource())
        assert 'must be ascii' in str(e.value)

    # make sure log gets warnings when we rename files
    upload(
        [x[0] for x in file_list], FakeSession(), FakeSamplesResource(), log=logging, coerce_ascii=True
    )

    for _, after in file_list:
        assert after in caplog.text


def test_single_end_files():
    session = FakeSession()
    samples_resource = FakeSamplesResource()
    upload(
        ['tests/data/files/test_R1_L001.fq.gz', 'tests/data/files/test_R2_L001.fq.gz'],
        session, samples_resource, threads=1
    )


def test_paired_end_files():
    session = FakeSession()
    samples_resource = FakeSamplesResource()
    upload(
        [('tests/data/files/test_R1_L001.fq.gz', 'tests/data/files/test_R2_L001.fq.gz')],
        session, samples_resource, threads=1
    )


def test_upload_fileobj():
    file_obj = BytesIO(b'>test\nACGT\n')
    session = FakeSession()
    samples_resource = FakeSamplesResource()
    upload_fileobj(file_obj, 'test.fa', 11, session, samples_resource)
    file_obj.close()


def test_multipart_encoder_passthru():
    wrapper = FASTXPassthru(
        'tests/data/files/test_R1_L001.fq.gz', os.path.getsize('tests/data/files/test_R1_L001.fq.gz')
    )
    wrapper_len = len(wrapper.read())
    wrapper.seek(0)

    assert wrapper_len == wrapper._fsize

    multipart_fields = OrderedDict()
    multipart_fields['file'] = ('fakefile', wrapper, 'application/x-gzip')
    encoder = MultipartEncoder(multipart_fields)
    MAGIC_HEADER_LEN = 178
    wrapper.seek(0)
    assert len(encoder.read()) - MAGIC_HEADER_LEN == wrapper_len


def test_multipart_encoder_interleave():
    pair = ('tests/data/files/test_R1_L001.fq.gz', 'tests/data/files/test_R2_L001.fq.gz')
    fname, fsize, fformat = _file_stats(pair)
    wrapper = FASTXInterleave(pair, fsize, fformat)
    wrappertext = wrapper.read()
    wrapper.seek(0)

    assert len(wrappertext) == fsize

    multipart_fields = OrderedDict()
    multipart_fields['file'] = ('fakefile', wrapper, 'text/plain')
    encoder = MultipartEncoder(multipart_fields)
    MAGIC_HEADER_LEN = 170  # shorter because of text/plain mime-type
    encodertext = encoder.read()
    assert len(encodertext) - MAGIC_HEADER_LEN == len(wrappertext)
