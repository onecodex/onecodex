"""
Unit tests for One Codex
"""
from pkg_resources import resource_filename
from onecodex.lib.sniff import sniff_file
from onecodex.lib.auth import check_version, fetch_api_key_from_uname
from onecodex.version import __version__


SERVER = 'https://app.onecodex.com/'


def test_sniffer():
    resp = sniff_file(resource_filename(__name__, 'data/files/test.fa'))

    assert resp['compression'] == 'none'
    assert resp['file_type'] == 'fasta'
    assert resp['seq_type'] == 'dna'

    assert not resp['seq_multiline']
    assert not resp['seq_has_gaps']
    assert not resp['seq_has_lowercase']
    assert not resp['seq_has_iupac']
    assert not resp['seq_has_unknowns']

    assert resp['seq_est_avg_len'] == 74.0
    assert resp['seq_est_gc'] == 0.5

    assert not resp['interleaved']

    resp = sniff_file(resource_filename(__name__, 'data/files/test.fq'))

    assert resp['compression'] == 'none'
    assert resp['file_type'] == 'fastq'
    assert resp['seq_type'] == 'dna'

    assert resp['qual_ids'] == 'blank_second'
    assert resp['qual_type'] == 'sanger'

    assert not resp['seq_multiline']
    assert not resp['seq_has_gaps']
    assert not resp['seq_has_lowercase']
    assert not resp['seq_has_iupac']
    assert resp['seq_has_unknowns']

    assert resp['seq_est_avg_len'] == 31.0

    assert not resp['interleaved']


def test_check_version():
    # TODO: Remove Internet dependency here -- need a version response mock
    should_upgrade, msg = check_version(__version__, SERVER, 'gui')

    assert not should_upgrade
    assert msg is None or msg.startswith('Please upgrade your client to the latest version')


def test_login():
    # TODO: Remove Internet dependency here -- need a version response mock
    #       May want a live tests suite that *does* do this though for awareness? Investigate.
    resp = fetch_api_key_from_uname('test', 'test', SERVER)
    assert resp is None  # should not be able to log in with this
