"""
Functions for connecting to the One Codex server; these should be used across CLI and GUI clients
"""
import re
import requests


class BearerTokenAuth(requests.auth.AuthBase):
    """Attaches Bearer Auth headers to a given Request object."""
    def __init__(self, token):
        self.token = token

    def __call__(self, request):
        request.headers['Authorization'] = "Bearer {}".format(self.token)
        return request


def fetch_api_key_from_uname(username, password, server_url):
    """
    Retrieves an API key from the One Codex webpage given the username and password
    """
    with requests.Session() as session:
        # get the login page normally
        text = session.get(server_url + 'login').text

        # retrieve the CSRF token out of it
        csrf = re.search('type="hidden" value="([^"]+)"', text).group(1)

        # and resubmit using the username/password *and* the CSRF
        login_data = {'email': username, 'password': password,
                      'csrf_token': csrf, 'next': '/api/get_token'}
        page = session.post(server_url + 'login', data=login_data)
        try:
            key = page.json()['key']
        except (ValueError, KeyError):  # ValueError includes simplejson.decoder.JSONDecodeError
            key = None
    return key


def check_version(version, server_url, client='cli'):
    """
    Check if the current version of the client software is supported by the One Codex
    backend. Returns a tuple with two values:
        - True if the user *must* upgrade their software, otherwise False
        - An error message if the user should upgrade, otherwise None.
    """
    def version_inadequate(client_version, server_version):
        """
        Simple, fast check for version inequality.

        Could use python package `semver` if we need more precise checks in
        edge cases, but this generally works for now.
        """
        client_version = tuple([int(x) for x in client_version.split('-')[0].split('.')])
        server_version = tuple([int(x) for x in server_version.split('.')])
        return client_version < server_version

    if client == 'cli':
        # TODO: update this for the new v1 client route
        data = requests.post(server_url + 'api/v0/check_for_cli_update', data={'version': version})
    elif client == 'gui':
        data = requests.post(server_url + 'api/v0/check_upload_app_version',
                             data={'version': version})
    else:
        raise Exception('Not a valid client descriptor')

    if data.status_code != 200:
        return False, 'Error connecting to server'
    data = data.json()
    latest_version = data['latest_version']

    if client == 'cli':
        uploader_text = ' using the command `pip install --upgrade onecodex`.'
    else:
        uploader_text = (' from the '
                         '<a href="http://www.onecodex.com/uploader.html">One Codex website</a>')

    # Only warn on minimum supported for version for the GUI, always warn about upgrades for CLI
    if client != 'cli' and 'min_supported_version' in data:
        min_version = data['min_supported_version']
        if version_inadequate(version, min_version):
            return True, ('Please upgrade your client to the latest version ' +
                          '(v{}){}; '.format(latest_version, uploader_text) +
                          'this version (v{}) is no longer supported.'.format(version))
        else:
            return False, None

    if version_inadequate(version, latest_version):
        return True, ('Please upgrade your client to the latest version ' +
                      '(v{}){}'.format(latest_version, uploader_text))

    return False, None
