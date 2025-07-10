"""
New API client implementation using httpx and Pydantic models.

This replaces the Potion-Client based implementation with a modern
HTTP client and auto-generated Pydantic models.
"""

import json
import logging
import os
import warnings
from typing import Dict, Optional

import filelock

from onecodex.client import OneCodexHTTPClient, ModelHTTPClient
from onecodex.exceptions import OneCodexException
from onecodex.lib.auth import BearerTokenAuth
from onecodex.models.base import OneCodexModel
from onecodex.utils import collapse_user, init_sentry
from onecodex.version import __version__

log = logging.getLogger("onecodex")


class Api:
    """One Codex API client using modern HTTP client and Pydantic models."""
    
    def __init__(
        self,
        api_key: Optional[str] = None,
        bearer_token: Optional[str] = None,
        base_url: Optional[str] = None,
        telemetry: Optional[bool] = None,
        load_extensions: bool = True,
        **kwargs
    ):
        """Initialize the One Codex API client.
        
        Args:
            api_key: API key for authentication
            bearer_token: Bearer token for authentication (preferred)
            base_url: Base URL for the API
            telemetry: Enable telemetry reporting
            load_extensions: Load visualization extensions
            **kwargs: Additional arguments
        """
        # Set default base URL
        if base_url is None:
            base_url = os.environ.get("ONE_CODEX_API_BASE", "https://app.onecodex.com")
            if base_url != "https://app.onecodex.com":
                warnings.warn(f"Using base API URL: {base_url}")
        
        self._base_url = base_url
        
        # Auto-fetch credentials if none provided
        if api_key is None and bearer_token is None:
            api_key, bearer_token = self._fetch_credentials()
        
        # Create HTTP client
        self._http_client = OneCodexHTTPClient(
            base_url=base_url,
            api_key=api_key,
            bearer_token=bearer_token,
            **kwargs
        )
        
        # Create model client
        self._model_client = ModelHTTPClient(self._http_client)
        
        # Bind models to this API client
        self._bind_models()
        
        # Configure visualization extensions
        if load_extensions:
            try:
                from onecodex.viz import configure_onecodex_theme
                configure_onecodex_theme()
            except ImportError:
                pass  # Visualization dependencies are optional
        
        # Configure telemetry
        if telemetry is True or (
            telemetry is None and os.environ.get("ONE_CODEX_AUTO_TELEMETRY", False)
        ):
            try:
                init_sentry(user_context={"email": self._fetch_account_email()})
                self._telemetry = True
            except Exception:
                self._telemetry = False
        else:
            self._telemetry = False
    
    def _fetch_credentials(self) -> tuple[Optional[str], Optional[str]]:
        """Fetch credentials from file or environment.
        
        Returns:
            Tuple of (api_key, bearer_token)
        """
        api_key = None
        bearer_token = None
        
        # Try to load from credentials file
        creds_file = os.path.expanduser("~/.onecodex")
        
        if os.path.exists(creds_file) and os.access(creds_file, os.R_OK):
            try:
                with filelock.FileLock(f"{creds_file}.lock"):
                    with open(creds_file, "r") as f:
                        creds = json.load(f)
                        api_key = creds.get("api_key")
            except (KeyError, ValueError, json.JSONDecodeError):
                warnings.warn(f"Credentials file ({collapse_user(creds_file)}) is corrupt")
        elif os.path.exists(creds_file):
            warnings.warn(f"Check permissions on {collapse_user(creds_file)}")
        
        # Try environment variables
        if api_key is None:
            api_key = os.environ.get("ONE_CODEX_API_KEY")
        if bearer_token is None:
            bearer_token = os.environ.get("ONE_CODEX_BEARER_TOKEN")
        
        return api_key, bearer_token
    
    def _fetch_account_email(self) -> Optional[str]:
        """Fetch account email for telemetry."""
        # Try credentials file first
        creds_file = os.path.expanduser("~/.onecodex")
        
        if os.path.exists(creds_file) and os.access(creds_file, os.R_OK):
            try:
                with open(creds_file, "r") as f:
                    creds = json.load(f)
                    return creds.get("email")
            except (KeyError, ValueError, json.JSONDecodeError):
                pass
        
        # Try environment variables
        return os.environ.get("ONE_CODEX_USER_EMAIL", os.environ.get("ONE_CODEX_USER_UUID"))
    
    def _bind_models(self):
        """Bind model classes to this API client."""
        from onecodex.models.samples import Samples, Metadata
        from onecodex.models.analyses import Analyses, Classifications, Alignments, FunctionalProfiles, Panels
        from onecodex.models.misc import Projects, Users, Jobs, Tags, Assets, Documents
        
        # Model registry
        model_classes = [
            Samples, Metadata, Analyses, Classifications, Alignments, 
            FunctionalProfiles, Panels, Projects, Users, Jobs, Tags, Assets, Documents
        ]
        
        # Bind each model class to this API client
        for model_class in model_classes:
            if hasattr(model_class, '_resource_path'):
                # Set class-level attributes
                model_class._api = self
                model_class._client = self._model_client
                
                # Add to this API instance
                setattr(self, model_class.__name__, model_class)
    
    def close(self):
        """Close the HTTP client."""
        self._http_client.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


# For backward compatibility, also export the new classes
# with the same names as the old ones
from onecodex.models.samples import Samples, Metadata
from onecodex.models.analyses import Analyses, Classifications, Alignments, FunctionalProfiles, Panels  
from onecodex.models.misc import Projects, Users, Jobs, Tags, Assets, Documents
from onecodex.models.new_collection import SampleCollection

__all__ = [
    "Api",
    "Samples", "Metadata", "Analyses", "Classifications", "Alignments",
    "FunctionalProfiles", "Panels", "Projects", "Users", "Jobs", "Tags", 
    "Assets", "Documents", "SampleCollection"
]