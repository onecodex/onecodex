from collections import OrderedDict
from io import BytesIO
from requests_toolbelt import MultipartEncoder

from mock import patch
import pytest

from onecodex.lib.inline_validator import FASTXTranslator
from onecodex.lib.upload import upload, upload_file, upload_large_file


@pytest.mark.parametrize('file_list,n_small,n_big', [
    (['file.1000.fa', 'file.5e10.fa'], 1, 1),
    ([('file.3e10.R1.fa', 'file.3e10.R2.fa')], 0, 1),
    ([('file.3e5.R1.fa', 'file.3e5.R2.fa'), 'file.1e5.fa'], 2, 0),
])
def test_upload_lots_of_files(file_list, n_small, n_big):
    session, samples_resource, server_url = None, None, None

    fake_size = lambda filename: int(float(filename.split('.')[1]))  # noqa

    uf = 'onecodex.lib.upload.upload_file'
    ulf = 'onecodex.lib.upload.upload_large_file'
    opg = 'onecodex.lib.upload.os.path.getsize'
    wf = 'onecodex.lib.upload._wrap_files'
    with patch(uf) as sm_upload, patch(ulf) as lg_upload, patch(wf) as p:
        with patch(opg, side_effect=fake_size) as p2:
            upload(file_list, session, samples_resource, server_url)
            assert sm_upload.call_count == n_small
            assert lg_upload.call_count == n_big
            assert p.call_count == len(file_list)
            assert p2.call_count == sum(2 if isinstance(f, tuple) else 1 for f in file_list)


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

    def read_init_multipart_upload(self):
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


def test_upload_big_file():
    file_obj = BytesIO(b'>test\nACGT\n')
    session = FakeSession()
    samples_resource = FakeSamplesResource()
    server_url = ''

    with patch('boto3.client') as p:
        upload_large_file(file_obj, 'test.fa', session,
                          samples_resource, server_url)
        assert p.call_count == 1

    file_obj.close()


def test_paired_end_upload():
    session = FakeSession()
    samples_resource = FakeSamplesResource()
    upload([('tests/data/files/test_R1_L001.fq.gz', 'tests/data/files/test_R2_L001.fq.gz')],
           session, samples_resource, '', threads=1)


def test_upload_small_file():
    file_obj = BytesIO(b'>test\nACGT\n')
    session = FakeSession()
    samples_resource = FakeSamplesResource()

    upload_file(file_obj, 'test.fa', session, samples_resource)
    file_obj.close()


def test_multipart_encoder():
    long_seq = b'ACGT' * 50
    data = b'\n'.join(b'>header_' + str(i).encode() + b'\n' + long_seq + b'\n'
                      for i in range(200))
    wrapper = FASTXTranslator(BytesIO(data), recompress=False)
    wrapper_len = len(wrapper.read())
    wrapper.seek(0)

    multipart_fields = OrderedDict()
    multipart_fields['file'] = ('fakefile', wrapper, 'application/x-gzip')
    encoder = MultipartEncoder(multipart_fields)
    MAGIC_HEADER_LEN = 178
    wrapper.seek(0)
    assert len(encoder.read()) - MAGIC_HEADER_LEN == wrapper_len
