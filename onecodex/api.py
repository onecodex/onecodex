from __future__ import print_function
import filelock
import json
import logging
import os
from requests.auth import HTTPBasicAuth
import warnings

from onecodex.lib.auth import BearerTokenAuth
from onecodex.utils import collapse_user, init_sentry, get_requests_session
from onecodex.version import __version__


log = logging.getLogger("onecodex")


class HTTPClient:
    """Simple requests-based HTTP client for API requests."""

    def __init__(self, auth=None, headers=None):
        self.session = get_requests_session(auth=auth, headers=headers)

    def get(self, url, **kwargs):
        """Make GET request."""
        return self.request("GET", url, **kwargs)

    def post(self, url, **kwargs):
        """Make POST request."""
        return self.request("POST", url, **kwargs)

    def patch(self, url, **kwargs):
        """Make PATCH request."""
        return self.request("PATCH", url, **kwargs)

    def put(self, url, **kwargs):
        """Make PUT request."""
        return self.request("PUT", url, **kwargs)

    def delete(self, url, **kwargs):
        """Make DELETE request."""
        return self.request("DELETE", url, **kwargs)

    def request(self, method, url, **kwargs):
        """Make request with specified method."""
        return self.session.request(method, url, **kwargs)


class Api(object):
    """One Codex Base API object class.

    Instantiates a Potion-Client object under the hood for making requests.
    """

    def __init__(
        self,
        api_key=None,
        bearer_token=None,
        base_url=None,
        telemetry=None,
        schema_path="/api/v1/schema",
        load_extensions=True,
        **kwargs,
    ):
        if base_url is None:
            base_url = os.environ.get("ONE_CODEX_API_BASE", "https://app.onecodex.com")
            if base_url != "https://app.onecodex.com":
                warnings.warn("Using base API URL: {}".format(base_url))

        self._base_url = base_url
        self._schema_path = schema_path

        # Attempt to automatically fetch API key from
        # ~/.onecodex file, API key, or bearer token environment vars
        # *if and only if* no auth is explicitly passed to Api
        #
        # TODO: Consider only doing this if an add'l env var like
        #       'ONE_CODEX_AUTO_LOGIN' or similar is set.
        if api_key is None and bearer_token is None:
            creds_file = os.path.expanduser("~/.onecodex")

            if not os.path.exists(creds_file):
                pass
            elif not os.access(creds_file, os.R_OK):
                warnings.warn("Check permissions on {}".format(collapse_user(creds_file)))
            else:
                try:
                    with filelock.FileLock("{}.lock".format(creds_file)):
                        api_key = json.load(open(creds_file, "r"))["api_key"]
                except KeyError:
                    # lacking an api_key doesn't mean the file is corrupt--it can just be that the
                    # schema was cached after logging in anonymously
                    pass
                except ValueError:
                    warnings.warn(
                        "Credentials file ({}) is corrupt".format(collapse_user(creds_file))
                    )

            if api_key is None:
                api_key = os.environ.get("ONE_CODEX_API_KEY")
            if bearer_token is None:
                bearer_token = os.environ.get("ONE_CODEX_BEARER_TOKEN")

        auth = None
        if bearer_token:  # prefer bearer token where available
            auth = BearerTokenAuth(bearer_token)
        elif api_key:
            auth = HTTPBasicAuth(api_key, "")

        headers = {"X-OneCodex-Client-User-Agent": __version__}

        if kwargs.get("experimental", False):
            warnings.warn(
                "Experimental API mode enabled. Features of the experimental API are subject to "
                "change without notice and should not be relied upon in a production environment."
            )
            headers["X-OneCodex-Api-Experimental"] = "1"

        self._client = HTTPClient(auth=auth, headers=headers)
        self._copy_resources()

        # Optionally configure custom One Codex altair theme and renderer
        if load_extensions:
            from onecodex.viz import configure_onecodex_theme

            configure_onecodex_theme()

        # Optionally configure Sentry
        if telemetry is True or (
            telemetry is None and os.environ.get("ONE_CODEX_AUTO_TELEMETRY", False)
        ):
            init_sentry(user_context={"email": self._fetch_account_email()})
            self._telemetry = True
        else:
            self._telemetry = False

    def _fetch_account_email(self):
        creds_file = os.path.expanduser("~/.onecodex")

        if not os.path.exists(creds_file):
            pass
        elif not os.access(creds_file, os.R_OK):
            warnings.warn("Check permissions on {}".format(collapse_user(creds_file)))
        else:
            try:
                return json.load(open(creds_file, "r"))["email"]
            except KeyError:
                pass
            except ValueError:
                warnings.warn("Credentials file ({}) is corrupt".format(collapse_user(creds_file)))

        return os.environ.get("ONE_CODEX_USER_EMAIL", os.environ.get("ONE_CODEX_USER_UUID"))

    def _copy_resources(self):
        """Copy models into this instance of Api().

        Returns
        -------
        `None`
        """
        from onecodex.models import _MODEL_REGISTRY

        for model in _MODEL_REGISTRY.values():
            model._api = self
            model._client = self._client
            assert getattr(self, model.__name__, None) is None
            setattr(self, model.__name__, model)
