from __future__ import print_function
import datetime
import json
import logging
import os
import platform
import requests
import sys
from requests.auth import HTTPBasicAuth
import warnings

from onecodex.lib.auth import (
    TOKEN_RENEWAL_WINDOW,
    BearerTokenAuth,
    InvalidCredentialsError,
    OneCodexCredentials,
)
from onecodex.utils import collapse_user, init_sentry, get_requests_session
from onecodex.version import __version__

from typing import TYPE_CHECKING, Type


log = logging.getLogger("onecodex")


def _read_creds_file() -> OneCodexCredentials | None:
    """Read the ~/.onecodex credentials file, returns None if it's unusable."""
    creds_file = os.path.expanduser("~/.onecodex")

    if not os.path.exists(creds_file):
        return None

    if not os.access(creds_file, os.R_OK):
        warnings.warn("Check permissions on {}".format(collapse_user(creds_file)))
        return None

    try:
        return OneCodexCredentials.load(creds_file)
    except InvalidCredentialsError:
        warnings.warn("Credentials file ({}) is corrupt".format(collapse_user(creds_file)))
        return None


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

    # Because models are added dynamically at runtime, static analysis tools will not know about
    # them. So we add some hints

    if TYPE_CHECKING:
        from onecodex.models import (
            Samples as _Samples,
            Metadata as _Metadata,
            Classifications as _Classifications,
            Analyses as _Analyses,
            Alignments as _Alignments,
            FunctionalProfiles as _FunctionalProfiles,
            Panels as _Panels,
            Mlsts as _Mlsts,
            Workflows as _Workflows,
            Users as _Users,
            Projects as _Projects,
            Tags as _Tags,
            Jobs as _Jobs,
            Documents as _Documents,
            Assets as _Assets,
            AnnotationSets as _AnnotationSets,
            Assemblies as _Assemblies,
            Genomes as _Genomes,
            Taxa as _Taxa,
        )

        Samples: Type[_Samples]
        Metadata: Type[_Metadata]
        Classifications: Type[_Classifications]
        Analyses: Type[_Analyses]
        Alignments: Type[_Alignments]
        FunctionalProfiles: Type[_FunctionalProfiles]
        Panels: Type[_Panels]
        Mlsts: Type[_Mlsts]
        Workflows: Type[_Workflows]
        Users: Type[_Users]
        Projects: Type[_Projects]
        Tags: Type[_Tags]
        Jobs: Type[_Jobs]
        Documents: Type[_Documents]
        Assets: Type[_Assets]
        AnnotationSets: Type[_AnnotationSets]
        Assemblies: Type[_Assemblies]
        Genomes: Type[_Genomes]
        Taxa: Type[_Taxa]

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

        # Attempt to automatically fetch auth from the bearer token or API key
        # environment vars, falling back to the login token or API key in the
        # ~/.onecodex file, *if and only if* no auth is explicitly passed to Api
        #
        # TODO: Consider only doing this if an add'l env var like
        #       'ONE_CODEX_AUTO_LOGIN' or similar is set.
        if api_key is None and bearer_token is None:
            bearer_token = os.environ.get("ONE_CODEX_BEARER_TOKEN")
            api_key = os.environ.get("ONE_CODEX_API_KEY")

            # The credentials file is used only without auth in the environment
            creds = None if bearer_token or api_key else _read_creds_file()
            if creds is not None:
                # Try bearer token from creds file first, fall back to api key.
                time_left = creds.time_left
                if time_left is None:
                    api_key = creds.api_key
                elif time_left <= datetime.timedelta(0):
                    warnings.warn(
                        "Your One Codex login has expired. Run `onecodex login` to log in again."
                    )
                else:
                    if time_left < TOKEN_RENEWAL_WINDOW:
                        warnings.warn(
                            "Your One Codex login expires soon. Run `onecodex login` to renew it."
                        )
                    bearer_token = creds.access_token

        auth = None
        if bearer_token:  # prefer bearer token where available
            auth = BearerTokenAuth(bearer_token)
        elif api_key:
            auth = HTTPBasicAuth(api_key, "")

        telemetry_enabled = telemetry is True or (
            telemetry is None and os.environ.get("ONE_CODEX_AUTO_TELEMETRY", False)
        )

        ocx_client_user_agent = {
            "onecodex_version": __version__,
            "requests_version": requests.__version__,
            "python_version": sys.version.split()[0] if telemetry_enabled else None,
            "os": platform.system() if telemetry_enabled else None,
            "os_version": platform.release() if telemetry_enabled else None,
        }

        headers = {
            "User-Agent": f"onecodex/{__version__}",
            "X-OneCodex-Client-User-Agent": json.dumps(ocx_client_user_agent),
        }

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
        if telemetry_enabled:
            init_sentry(user_context={"email": self._fetch_account_email()})
            self._telemetry = True
        else:
            self._telemetry = False

    def _fetch_account_email(self):
        creds = _read_creds_file()
        if creds is not None and creds.email:
            return creds.email

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
