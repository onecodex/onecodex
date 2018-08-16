"""
Unit tests for One Codex
"""
from onecodex.lib.auth import check_version
from onecodex.version import __version__


SERVER = 'https://app.onecodex.com/'


def test_check_version_integration():
    # TODO: Remove Internet dependency here -- need a version response mock
    should_upgrade, msg = check_version(__version__, SERVER, 'cli')

    assert not should_upgrade
    assert msg is None or msg.startswith('Please upgrade your client to the latest version')
