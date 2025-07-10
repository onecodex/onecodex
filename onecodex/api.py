from __future__ import print_function
import filelock
import json
import logging
import os
import requests
from requests.adapters import HTTPAdapter
from requests.auth import HTTPBasicAuth
from requests.packages.urllib3.util.retry import Retry
import warnings

from onecodex.lib.auth import BearerTokenAuth
from onecodex.utils import collapse_user, init_sentry
from onecodex.version import __version__


log = logging.getLogger("onecodex")


class HTTPClient:
    """Simple requests-based HTTP client for API requests."""

    def __init__(self, auth=None, headers=None, timeout=30.0, **kwargs):
        self.session = requests.Session()
        if auth:
            self.session.auth = auth
        if headers:
            self.session.headers.update(headers)

        # Set backoff / retry strategy for 429s
        retry_strategy = Retry(
            total=3, backoff_factor=4, allowed_methods=None, status_forcelist=[429, 502, 503]
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)

    def get(self, url, **kwargs):
        """Make GET request."""
        return self.session.get(url, **kwargs)

    def post(self, url, **kwargs):
        """Make POST request."""
        return self.session.post(url, **kwargs)

    def put(self, url, **kwargs):
        """Make PUT request."""
        return self.session.put(url, **kwargs)

    def delete(self, url, **kwargs):
        """Make DELETE request."""
        return self.session.delete(url, **kwargs)

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
        cache_schema=True,
        base_url=None,
        telemetry=None,
        schema_path="/api/v1/schema",
        load_extensions=True,
        **kwargs,
    ):
        if cache_schema is False:
            warnings.warn(
                "The cache_schema parameter is deprecated and will be removed in a future version.",
                DeprecationWarning,
                stacklevel=2,
            )
        if base_url is None:
            base_url = os.environ.get("ONE_CODEX_API_BASE", "https://app.onecodex.com")
            if base_url != "https://app.onecodex.com":
                warnings.warn("Using base API URL: {}".format(base_url))

        if kwargs.get("experimental", False):
            warnings.warn(
                "Experimental API mode enabled. Features of the experimental API are subject to "
                "change without notice and should not be relied upon in a production environment."
            )
            schema_path = "/api/v1_experimental/schema"
            cache_schema = False

        self._req_args = {}
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

        self._http_client = HTTPClient(auth=auth, headers=headers)
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
        """Copy subclassed potion resources into this instance of Api().

        Notes
        -----
        Will only make available any objects that were listed in the REST API schema, either cached
        or retrieved from the server. If a non-standard schema_path was given, make sure we match
        e.g., /api/v1_experimental/samples up with /api/v1/samples.

        Returns
        -------
        `None`
        """
        pass

    #     from onecodex.models import _model_lookup

    #     if self._schema_path != "/api/v1/schema":
    #         new_schema_base = self._schema_path.rstrip("schema")

    #         aliased_resources = {}

    #         for resource, oc_cls in _model_lookup.items():
    #             aliased_resources[resource.replace("/api/v1/", new_schema_base)] = oc_cls

    #         _model_lookup.update(aliased_resources)

    #     for resource in self._client._resources:
    #         # set the name param, the keys now have / in them
    #         potion_resource = self._client._resources[resource]

    #         try:
    #             oc_cls = _model_lookup[resource]
    #             oc_cls._api = self
    #             oc_cls._resource = potion_resource
    #             setattr(self, oc_cls.__name__, oc_cls)
    #         except KeyError:  # Ignore resources we don't explicitly model
    #             pass
