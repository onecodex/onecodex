from __future__ import print_function
import filelock
import json
import logging
import os
import platform
import requests
import sys
from requests.auth import HTTPBasicAuth
import warnings

from onecodex.lib.auth import BearerTokenAuth
from onecodex.utils import collapse_user, init_sentry, get_requests_session
from onecodex.version import __version__

from typing import TYPE_CHECKING, Type


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
    """Entry point to the One Codex API.

    Instantiating `Api()` wires up authentication, attaches every model class
    (`Samples`, `Classifications`, `FunctionalProfiles`, ...) to this client, and
    optionally configures a disk cache for analysis results.

    Authentication is resolved in this order, stopping at the first match:
        1. `api_key` / `bearer_token` passed explicitly.
        2. `~/.onecodex` credentials file (written by `onecodex login`).
        3. `ONE_CODEX_API_KEY` / `ONE_CODEX_BEARER_TOKEN` env vars.

    Most users want the default::

        from onecodex import Api
        ocx = Api()
        sample = ocx.Samples.get("...")

    Parameters
    ----------
    api_key : str, optional
        API key for HTTP basic auth. Falls back to creds file / env var.
    bearer_token : str, optional
        Bearer token. Preferred over `api_key` when both are available.
    base_url : str, optional
        API base URL. Defaults to `https://app.onecodex.com`, overridable via
        the `ONE_CODEX_API_BASE` env var.
    telemetry : bool, optional
        Enable Sentry telemetry. Defaults to off unless `ONE_CODEX_AUTO_TELEMETRY`
        is set.
    schema_path : str, optional
        Path to the API schema endpoint.
    load_extensions : bool, optional
        If True (default), install the One Codex Altair theme/renderer.
    cache_results : bool or str, optional
        Opt-in disk cache for analysis `.results()` / `.table()` fetches across
        processes — useful when the same analysis is fetched repeatedly across sessions.
        Defaults to False (no cache).

        - `False` (default): no caching.
        - `True`: cache under the system tempdir (per-user, ephemeral — the OS
          reclaims it on its own schedule, so disk fill isn't a concern).
        - `str`: an explicit path. Use this when you want the cache to persist
          across reboots / sandbox runs. You manage the disk budget.

        The `ONECODEX_DISK_CACHE` env var sets a path when `cache_results` is
        left at the default. See `onecodex.cache.DiskCache` for the on-disk
        format and freshness semantics.
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
        api_key: str | None = None,
        bearer_token: str | None = None,
        base_url: str | None = None,
        telemetry: None | bool = None,
        schema_path: str = "/api/v1/schema",
        load_extensions: bool = True,
        cache_results: bool | str = False,
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

        env_cache = os.environ.get("ONECODEX_DISK_CACHE")
        if cache_results is False and env_cache:
            cache_results = env_cache
        if cache_results:
            from onecodex.cache import DiskCache

            self._cache = DiskCache(path=cache_results if isinstance(cache_results, str) else None)
        else:
            self._cache = None

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
