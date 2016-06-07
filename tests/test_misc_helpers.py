from onecodex.lib.auth import check_version
from onecodex.lib.sniff import sniff
import pytest
from tests.conftest import mock_requests


@pytest.mark.parametrize('client_version,server_version,min_version,client,requires_update', [
    ('0.1.3', '0.1.3', None, 'cli', False),
    ('0.1.4', '0.1.3', None, 'cli', False),
    ('0.1.3', '0.1.4', None, 'cli', True),
    ('0.2.3', '0.10.4', None, 'cli', True),
    ('0.1.1', '0.1.3', '0.1.2', 'gui', True),
    ('0.1.2', '0.1.3', '0.1.2', 'gui', False),
    ('0.1.1', '0.2.0', '0.2.0', 'gui', True),
    ('0.2.0', '0.2.0', '0.2.0', 'gui', False),
])
def test_check_version(client_version, server_version, min_version, client, requires_update):
    json_data = {
        "POST::api/v0/check_for_cli_update": {
            "latest_version": server_version
        },
        "POST::api/v0/check_upload_app_version": {
            "latest_version": server_version,
            "min_supported_version": min_version
        }
    }
    url = "http://localhost:3000/"
    with mock_requests(json_data):
        result, error_code = check_version(client_version, url, client=client)
        assert result == requires_update


@pytest.mark.parametrize(('data,file_type,seq_multiline,seq_has_lowercase,'
                          'seq_est_gc,interleaved,qual_type'), [
    ('>Header\nAGAGAGAGAGAGAGAG\n', 'fasta', False, False, 0.5, False, None),
    ('>Header\nAGAGAGAGAGAGAGAG\nAGAGAGAGAGAGAGAG\n', 'fasta', True, False, 0.5, False, None),
    ('>Header\nAGAGAGAGAGAGAGAG\nAGAGAGAGAgaGAGAG\n', 'fasta', True, True, 0.5, False, None),
    # Single line FASTQ fails
    ('@Header1\nAGAGAGAGAGAGAGAG\n+\naaaaaaaaaaaaaaaa\n', 'bad',
     False, False, 0.5, False, 'sanger'),
    # Two line FASTQ passes
    ('@Header1\nAGAG\n+\naaaa\n@Header2\nGAGA\n+\naaaa\n', 'fastq',
     False, False, 0.5, False, 'sanger'),
    # Bad qual score, note we only parse the first one!?
    ('@Header1\nAGAG\n+\n!x!x\n@Header2\nGAGA\n+\n!x!x\n', 'fastq',
     False, False, 0.5, False, 'bad'),
])
def test_sniffer(data, file_type, seq_multiline, seq_has_lowercase,
                 seq_est_gc, interleaved, qual_type):
    start = data[0]
    data = data[1:]
    fastx_info = sniff(start, data)
    print(data)
    print(fastx_info)
    assert fastx_info['file_type'] == file_type
    if file_type != 'bad':
        assert fastx_info['seq_multiline'] == seq_multiline
        assert fastx_info['seq_has_lowercase'] == seq_has_lowercase
        assert fastx_info['seq_est_gc'] == seq_est_gc
        assert fastx_info['interleaved'] == interleaved
