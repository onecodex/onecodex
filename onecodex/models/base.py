"""
Base model classes that extend the auto-generated Pydantic models.

This module provides the foundation for One Codex models by extending
the auto-generated models with One Codex-specific functionality.
"""

import json
import logging
from datetime import datetime
from typing import Any, Dict, List, Optional, Type, TypeVar, Union

import pytz
from dateutil.parser import parse
from pydantic import BaseModel, Field, model_validator

from onecodex.client import ModelHTTPClient, OneCodexHTTPClient
from onecodex.exceptions import MethodNotSupported, OneCodexException, ServerError

log = logging.getLogger("onecodex")

ModelType = TypeVar("ModelType", bound="OneCodexModel")


class OneCodexModel(BaseModel):
    """Base class for all One Codex models.

    Provides common functionality like CRUD operations, property access,
    and API interactions that were previously handled by OneCodexBase.
    """

    # Class-level attributes (set by API client)
    _api: Optional["Api"] = None  # Forward reference to avoid circular imports
    _client: Optional[ModelHTTPClient] = None
    _resource_path: str = None  # Default resource path, subclasses override

    # Instance attributes
    id: Optional[str] = Field(None, alias="$uri")  # Map $uri to id for compatibility

    model_config = {
        "extra": "allow",  # Allow extra fields from API
        "populate_by_name": True,  # Allow field aliases
        "json_encoders": {datetime: lambda v: v.isoformat() if v else None},
    }

    def __init__(self, **data):
        """Initialize model with data conversion."""
        # Convert datetime strings to datetime objects
        converted_data = self._convert_datetime_fields(data)
        super().__init__(**converted_data)

    def __repr__(self) -> str:
        """String representation of the model."""
        return f"<{self.__class__.__name__} {self.id}>"

    def _convert_datetime_fields(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Convert datetime string fields to datetime objects."""
        converted = data.copy()

        # Get model fields that should be datetime
        datetime_fields = []
        for field_name, field_info in self.model_fields.items():
            if field_info.annotation == datetime:
                datetime_fields.append(field_name)

        for field_name in datetime_fields:
            if field_name in converted and isinstance(converted[field_name], str):
                try:
                    dt = parse(converted[field_name])
                    if dt.tzinfo is None:
                        dt = pytz.utc.localize(dt)
                    else:
                        dt = dt.astimezone(pytz.utc)
                    converted[field_name] = dt
                except (ValueError, TypeError):
                    pass  # Keep original value if parsing fails

        return converted

    def to_json(self, include_references: bool = True) -> str:
        """Convert model to JSON string.

        Args:
            include_references: Whether to include reference fields

        Returns:
            JSON string representation
        """
        data = self.model_dump(by_alias=True)

        if not include_references:
            # Remove reference fields (those starting with $ or containing $ref)
            data = {k: v for k, v in data.items() if not k.startswith("$") and "$ref" not in str(v)}

        return json.dumps(data, default=str)

    @classmethod
    def _check_bound(cls):
        """Check if the model is bound to an API client."""
        if cls._api is None or cls._client is None:
            raise OneCodexException(f"{cls.__name__} is not bound to an API client")

    @classmethod
    def _convert_id_to_uri(cls, uuid: str) -> str:
        """Convert UUID to full URI format.

        Args:
            uuid: Resource UUID

        Returns:
            Full URI for the resource
        """
        if not cls._resource_path:
            raise OneCodexException(f"No resource path defined for {cls.__name__}")

        if uuid.startswith(cls._resource_path):
            return uuid
        return f"{cls._resource_path}/{uuid}"

    @classmethod
    def get(cls: Type[ModelType], uuid: str) -> Optional[ModelType]:
        """Retrieve a model by UUID.

        Args:
            uuid: Resource UUID

        Returns:
            Model instance or None if not found
        """
        cls._check_bound()

        try:
            return cls._client.get_model(cls._resource_path, cls, uuid)
        except Exception as e:
            log.debug(f"Error fetching {cls.__name__} {uuid}: {e}")
            return None

    @classmethod
    def where(
        cls: Type[ModelType],
        sort: Optional[Union[str, List[str]]] = None,
        limit: Optional[int] = None,
        **filters,
    ) -> List[ModelType]:
        """Filter and retrieve models.

        Args:
            sort: Sort field(s)
            limit: Maximum number of results
            **filters: Filter parameters

        Returns:
            List of matching models
        """
        cls._check_bound()

        params = {}

        # Add filters
        if filters:
            # Remove special parameters
            filter_func = filters.pop("filter", None)
            public = filters.pop("public", False)

            # Add remaining filters to params
            params.update(filters)

        # Add sort parameter
        if sort:
            if isinstance(sort, list):
                params["sort"] = ",".join(sort)
            else:
                params["sort"] = sort

        # Add limit
        if limit:
            params["limit"] = limit

        # Make request
        models = cls._client.list_models(cls._resource_path, cls, params)

        # Apply local filter function if provided
        if "filter" in filters and callable(filters["filter"]):
            models = [m for m in models if filters["filter"](m)]

        return models

    @classmethod
    def all(
        cls: Type[ModelType],
        sort: Optional[Union[str, List[str]]] = None,
        limit: Optional[int] = None,
    ) -> List[ModelType]:
        """Retrieve all models of this type.

        Args:
            sort: Sort field(s)
            limit: Maximum number of results

        Returns:
            List of all models
        """
        return cls.where(sort=sort, limit=limit)

    def save(self) -> "OneCodexModel":
        """Save the model (create or update).

        Returns:
            Updated model instance
        """
        self._check_bound()

        if self.id is None:
            # Create new resource
            data = self.model_dump(by_alias=True, exclude_unset=True)
            return self._client.create_model(self._resource_path, self.__class__, data)
        else:
            # Update existing resource
            data = self.model_dump(by_alias=True, exclude_unset=True)
            return self._client.update_model(self._resource_path, self.__class__, self.id, data)

    def delete(self) -> bool:
        """Delete the model.

        Returns:
            True if deletion was successful
        """
        self._check_bound()

        if self.id is None:
            raise ServerError(f"{self.__class__.__name__} object does not exist yet")

        return self._client.delete_model(self._resource_path, self.id)


class ResourceList:
    """List-like container for One Codex model objects.

    Provides compatibility with the existing ResourceList functionality
    while working with Pydantic models.
    """

    def __init__(self, items: List[OneCodexModel], model_class: Type[OneCodexModel]):
        """Initialize with a list of model instances.

        Args:
            items: List of model instances
            model_class: The model class for type checking
        """
        self._items = items
        self._model_class = model_class

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self):
        return iter(self._items)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return ResourceList(self._items[index], self._model_class)
        return self._items[index]

    def __setitem__(self, index, value):
        if not isinstance(value, self._model_class):
            raise ValueError(f"Expected {self._model_class.__name__}, got {type(value).__name__}")
        self._items[index] = value

    def __delitem__(self, index):
        del self._items[index]

    def __contains__(self, item) -> bool:
        return item in self._items

    def __repr__(self) -> str:
        return f"ResourceList({self._items})"

    def append(self, item: OneCodexModel):
        """Add an item to the list."""
        if not isinstance(item, self._model_class):
            raise ValueError(f"Expected {self._model_class.__name__}, got {type(item).__name__}")
        self._items.append(item)

    def extend(self, items: List[OneCodexModel]):
        """Extend the list with multiple items."""
        for item in items:
            if not isinstance(item, self._model_class):
                raise ValueError(
                    f"Expected {self._model_class.__name__}, got {type(item).__name__}"
                )
        self._items.extend(items)

    def remove(self, item: OneCodexModel):
        """Remove an item from the list."""
        self._items.remove(item)

    def copy(self) -> "ResourceList":
        """Create a copy of the list."""
        return ResourceList(self._items.copy(), self._model_class)
