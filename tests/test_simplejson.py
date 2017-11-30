import pytest  # noqa

from onecodex import Cli


def test_simplejson(runner, upload_mocks):
    result = runner.invoke(Cli, ['upload'])
    assert result.exit_code != 0
    assert "You currently have simplejson installed" in result.output
