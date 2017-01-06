from onecodex.lib.auth import check_version
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
