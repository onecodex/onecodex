"""
HTTP client using requests library for backward compatibility with tests.
"""

import json
import logging
import os
from typing import Any, Dict, List, Optional, Type, TypeVar, Union
from urllib.parse import urljoin

import requests
from requests.adapters import HTTPAdapter
from requests.auth import HTTPBasicAuth
from requests.packages.urllib3.util.retry import Retry
from pydantic import BaseModel

from onecodex.exceptions import OneCodexException, PermissionDenied, ServerError
from onecodex.lib.auth import BearerTokenAuth
from onecodex.version import __version__

log = logging.getLogger("onecodex")

ModelType = TypeVar("ModelType", bound=BaseModel)


class OneCodexHTTPClient:
    """HTTP client for One Codex API using requests library."""

    def __init__(
        self,
        base_url: str = "https://app.onecodex.com",
        api_key: Optional[str] = None,
        bearer_token: Optional[str] = None,
        timeout: float = 30.0,
        retries: int = 3,
        **kwargs,
    ):
        """Initialize the HTTP client.

        Args:
            base_url: Base URL for the One Codex API
            api_key: API key for authentication
            bearer_token: Bearer token for authentication (preferred over api_key)
            timeout: Request timeout in seconds
            retries: Number of retries for failed requests
        """
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.retries = retries

        # Create session
        self.session = requests.Session()

        # Set up authentication
        if bearer_token:
            self.session.auth = BearerTokenAuth(bearer_token)
        elif api_key:
            self.session.auth = HTTPBasicAuth(api_key, "")

        # Set up headers
        self.session.headers.update(
            {
                "X-OneCodex-Client-User-Agent": __version__,
                "Content-Type": "application/json",
                "Accept": "application/json",
            }
        )

        # Configure retry strategy
        retry_strategy = Retry(
            total=retries, backoff_factor=4, allowed_methods=None, status_forcelist=[429, 502, 503]
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        """Close the HTTP client."""
        self.session.close()

    def _handle_response(self, response: requests.Response) -> Optional[Dict[str, Any]]:
        """Handle HTTP response and convert to JSON.

        Args:
            response: requests Response object

        Returns:
            Parsed JSON response or None for 404

        Raises:
            OneCodexException: For various API errors
        """
        if response.status_code == 404:
            return None

        if response.status_code == 401:
            raise PermissionDenied("Authentication failed")
        elif response.status_code == 403:
            raise PermissionDenied("Access forbidden")
        elif response.status_code == 400:
            try:
                error_data = response.json()
                error_msg = self._format_error_message(error_data)
                raise ServerError(error_msg)
            except (json.JSONDecodeError, KeyError):
                raise ServerError("Bad request")
        elif response.status_code == 409:
            raise ServerError("Resource conflict")
        elif response.status_code >= 500:
            raise ServerError(f"Server error: {response.status_code}")
        elif response.status_code >= 400:
            raise OneCodexException(f"HTTP {response.status_code}: {response.text}")

        try:
            return response.json()
        except json.JSONDecodeError:
            if response.status_code == 204:  # No Content
                return {}
            raise OneCodexException("Invalid JSON response")

    def _format_error_message(self, error_data: Dict[str, Any]) -> str:
        """Format error message from API response."""
        if isinstance(error_data, dict):
            errors = error_data.get("errors", [])
            if errors and isinstance(errors, list):
                # Handle validation errors
                if len(errors) == 1 and "validationOf" in errors[0]:
                    required = errors[0]["validationOf"].get("required", [])
                    return f"Validation error. Required fields: {', '.join(required)}"

                # Handle general errors
                messages = [err.get("message", str(err)) for err in errors]
                return "; ".join(filter(None, messages))

            # Single error message
            if "message" in error_data:
                return error_data["message"]

        return "Bad request"

    def get(
        self, endpoint: str, params: Optional[Dict[str, Any]] = None
    ) -> Optional[Dict[str, Any]]:
        """Make a GET request.

        Args:
            endpoint: API endpoint (e.g., "/api/v1/samples")
            params: Query parameters

        Returns:
            JSON response or None if 404
        """
        url = urljoin(self.base_url + "/", endpoint.lstrip("/"))
        response = self.session.get(url, params=params, timeout=self.timeout)
        return self._handle_response(response)

    def post(self, endpoint: str, data: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """Make a POST request.

        Args:
            endpoint: API endpoint
            data: Request body data

        Returns:
            JSON response
        """
        url = urljoin(self.base_url + "/", endpoint.lstrip("/"))
        response = self.session.post(url, json=data, timeout=self.timeout)
        return self._handle_response(response)

    def put(self, endpoint: str, data: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """Make a PUT request.

        Args:
            endpoint: API endpoint
            data: Request body data

        Returns:
            JSON response
        """
        url = urljoin(self.base_url + "/", endpoint.lstrip("/"))
        response = self.session.put(url, json=data, timeout=self.timeout)
        return self._handle_response(response)

    def patch(self, endpoint: str, data: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """Make a PATCH request.

        Args:
            endpoint: API endpoint
            data: Request body data

        Returns:
            JSON response
        """
        url = urljoin(self.base_url + "/", endpoint.lstrip("/"))
        response = self.session.patch(url, json=data, timeout=self.timeout)
        return self._handle_response(response)

    def delete(self, endpoint: str) -> Optional[Dict[str, Any]]:
        """Make a DELETE request.

        Args:
            endpoint: API endpoint

        Returns:
            JSON response or None
        """
        url = urljoin(self.base_url + "/", endpoint.lstrip("/"))
        response = self.session.delete(url, timeout=self.timeout)
        return self._handle_response(response)

    def paginated_get(
        self, endpoint: str, params: Optional[Dict[str, Any]] = None, page_size: int = 200
    ) -> List[Dict[str, Any]]:
        """Make paginated GET requests to fetch all results.

        Args:
            endpoint: API endpoint
            params: Query parameters
            page_size: Number of items per page

        Returns:
            List of all items from all pages
        """
        all_items = []
        current_params = (params or {}).copy()
        current_params["per_page"] = page_size
        offset = 0

        while True:
            current_params["offset"] = offset
            response = self.get(endpoint, current_params)

            if not response:
                break

            # Handle different response formats
            items = response.get("results", response) if isinstance(response, dict) else response

            if not items or not isinstance(items, list):
                break

            all_items.extend(items)

            # Check if we got fewer items than requested (end of results)
            if len(items) < page_size:
                break

            offset += page_size

        return all_items


class ModelHTTPClient:
    """HTTP client wrapper that works with Pydantic models."""

    def __init__(self, http_client: OneCodexHTTPClient):
        """Initialize with an HTTP client.

        Args:
            http_client: Configured OneCodexHTTPClient instance
        """
        self.http_client = http_client

    def get_model(
        self, endpoint: str, model_class: Type[ModelType], resource_id: str
    ) -> Optional[ModelType]:
        """Fetch a single model by ID.

        Args:
            endpoint: Base API endpoint (e.g., "/api/v1/samples")
            model_class: Pydantic model class
            resource_id: Resource ID

        Returns:
            Model instance or None if not found
        """
        url = f"{endpoint.rstrip('/')}/{resource_id}"
        data = self.http_client.get(url)

        if data is None:
            return None

        return model_class.model_validate(data)

    def list_models(
        self, endpoint: str, model_class: Type[ModelType], params: Optional[Dict[str, Any]] = None
    ) -> List[ModelType]:
        """Fetch a list of models.

        Args:
            endpoint: API endpoint
            model_class: Pydantic model class
            params: Query parameters

        Returns:
            List of model instances
        """
        items = self.http_client.paginated_get(endpoint, params)
        return [model_class.model_validate(item) for item in items]

    def create_model(
        self, endpoint: str, model_class: Type[ModelType], data: Dict[str, Any]
    ) -> ModelType:
        """Create a new model.

        Args:
            endpoint: API endpoint
            model_class: Pydantic model class
            data: Creation data

        Returns:
            Created model instance
        """
        response_data = self.http_client.post(endpoint, data)
        return model_class.model_validate(response_data)

    def update_model(
        self, endpoint: str, model_class: Type[ModelType], resource_id: str, data: Dict[str, Any]
    ) -> ModelType:
        """Update an existing model.

        Args:
            endpoint: Base API endpoint
            model_class: Pydantic model class
            resource_id: Resource ID
            data: Update data

        Returns:
            Updated model instance
        """
        url = f"{endpoint.rstrip('/')}/{resource_id}"
        response_data = self.http_client.patch(url, data)
        return model_class.model_validate(response_data)

    def delete_model(self, endpoint: str, resource_id: str) -> bool:
        """Delete a model.

        Args:
            endpoint: Base API endpoint
            resource_id: Resource ID

        Returns:
            True if deleted successfully
        """
        url = f"{endpoint.rstrip('/')}/{resource_id}"
        self.http_client.delete(url)
        return True
