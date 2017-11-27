import pytest

from click.testing import CliRunner

from onecodex import Cli

def test_simplejson(runner):
    result = runner.invoke(Cli, ['upload'])
    assert result.exit_code == 1
    assert "You currently have simplejson installed" in result.output
