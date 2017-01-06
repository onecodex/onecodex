"""
Unit tests for One Codex
"""
from onecodex.lib.auth import check_version, fetch_api_key_from_uname
from onecodex.version import __version__


SERVER = 'https://app.onecodex.com/'


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
