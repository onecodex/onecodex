import re
import requests


class BearerTokenAuth(requests.auth.AuthBase):
    """Attache Bearer Auth headers to the given Request object."""

    def __init__(self, token):
        self.token = token

    def __call__(self, request):
        request.headers["Authorization"] = "Bearer {}".format(self.token)
        return request


def fetch_api_key_from_uname(username, password, server_url):
    """Retrieve an API key from the One Codex webpage given the username and password."""
    # TODO: Hit programmatic endpoint to fetch JWT key, not API key
    # TODO: Open browser for a login and redirect to the client.
    #       This would enable a workflow with SSO.
    with requests.Session() as session:
        # get the login page normally
        text = session.get(server_url + "login").text

        # retrieve the CSRF token out of it
        csrf = re.search('type="hidden" value="([^"]+)"', text).group(1)

        # and resubmit using the username/password *and* the CSRF
        login_data = {
            "email": username,
            "password": password,
            "csrf_token": csrf,
            "next": "/api/get_token",
        }
        page = session.post(server_url + "login", data=login_data)
        try:
            key = page.json()["key"]
        except (ValueError, KeyError):  # ValueError includes simplejson.decoder.JSONDecodeError
            key = None
    return key


def check_version(version, server):
    """Check if the current CLI version is supported by the One Codex backend.

    Parameters
    ----------
    version : `string`
        Current version of the One Codex client library
    server : `string`
        Complete URL to One Codex server, e.g., https://app.onecodex.com

    Returns
    -------
    `tuple` containing two values:
        - True if the user *must* upgrade their software, otherwise False
        - An error message if the user should upgrade, otherwise None.
    """

    def version_inadequate(client_version, server_version):
        """Check for version inequality.

        Could use python package `semver` if we need more precise checks in edge cases, but this
        generally works for now.
        """
        client_version = tuple([int(x) for x in client_version.split("-")[0].split(".")])
        server_version = tuple([int(x) for x in server_version.split(".")])
        return client_version < server_version

    # this will probably live on /api/v0 forever for compat with older CLI versions
    data = requests.post(server + "api/v0/check_for_cli_update", data={"version": version})

    if data.status_code != 200:
        return False, "Error connecting to server"

    data = data.json()
    latest_version = data["latest_version"]

    if version_inadequate(version, latest_version):
        return (
            True,
            (
                "Please upgrade your client to the latest version (v{}) using the command "
                "`pip install --upgrade onecodex`".format(latest_version)
            ),
        )

    return False, None
