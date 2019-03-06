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
    # TODO: Hit programmatic endpoint to fetch JWT key, not API key
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
