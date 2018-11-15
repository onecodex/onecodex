from requests.auth import AuthBase


class HTTPBearerAuth(AuthBase):
    """Attaches HTTP Basic Authentication to the given Request object."""
    def __init__(self, token):
        self.token = token

    def __call__(self, r):
        r.headers['Authorization'] = 'Bearer {}'.format(self.token)
        return r
