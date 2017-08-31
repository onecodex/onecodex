import os

import mock
import pytest
from raven import Client as RavenClient

from onecodex import Cli
from onecodex.utils import get_raven_client, telemetry


def test_raven_off_by_default(ocx):
    assert ocx._telemetry is False
    assert ocx._raven_client is None


def test_no_raven_in_test_env(runner):
    assert runner.env.get('ONE_CODEX_NO_TELEMETRY') is not None


@pytest.mark.parametrize('cli_args,call_count', [
    (['logout'], 2),  # Once for the `onecodex` entry and once for `logout` since no Api client
    (['--no-telemetry', 'logout'], 0),
    (['--version'], 0)  # Built-in, no @click.pass_context
])
def test_telemetry_decorator(cli_args, call_count):

    @telemetry
    def cli_exception(cli_args):
        Cli.main(cli_args)

    with mock.patch('onecodex.utils.get_raven_client') as grc:
        with pytest.raises(SystemExit):
            cli_exception(cli_args)
        assert grc.call_count == call_count


def test_ocx_with_raven(ocx_w_raven, mocked_creds_file):
    assert ocx_w_raven._telemetry is True
    assert isinstance(ocx_w_raven._raven_client, RavenClient)
    assert ocx_w_raven._raven_client.remote.base_url == 'https://sentry.example.com'
    assert ocx_w_raven._raven_client.is_enabled() is True
    assert 'email' in ocx_w_raven._raven_client.context['user']


def test_good_sentry_dsn():
    raven = get_raven_client()
    assert isinstance(raven, RavenClient)


def test_bad_sentry_dsn():
    with mock.patch.dict(os.environ, {'ONE_CODEX_SENTRY_DSN': 'bad'}):
        raven = get_raven_client()
        assert raven is None
