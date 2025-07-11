from typing import ClassVar, Optional, List

from pydantic import BaseModel as PydanticBaseModel
from pydantic import ConfigDict, computed_field, Field


class ApiRef(PydanticBaseModel):
    """Represents an unresolved API reference."""

    ref: str = Field(..., alias="$ref")

    @computed_field
    @property
    def id(self) -> str:
        """Extract the ID from the reference URI."""
        return self.ref.split("/")[-1]

    @computed_field
    @property
    def resource_type(self) -> str:
        """Extract the resource type from the reference URI."""
        parts = self.ref.split("/")
        if len(parts) >= 4 and parts[1] == "api" and parts[2] == "v1":
            return parts[3]
        return "unknown"

    def __repr__(self):
        return f"<ApiRef {self.resource_type}:{self.id}>"

    def _resolve(self):
        """Resolve this reference to an actual object."""
        # Make API call to resolve
        target_class = self._get_target_class()
        if not target_class:
            raise ValueError(f"Cannot determine target class for {self.ref}")

        resp = target_class._client.get(
            f"{target_class._api._base_url}{target_class._resource_path}/{self.id}"
        )
        return target_class.model_validate(resp.json())

    def _get_target_class(self):
        """Determine the target class from the reference URI."""
        from onecodex.models.pydantic import get_model_class

        return get_model_class(self.ref)


class ApiBaseModel(PydanticBaseModel):
    _api: ClassVar[Optional["Api"]] = (  # noqa: F821
        None  # Forward reference to avoid circular imports
    )
    _client: ClassVar[Optional["HTTPClient"]] = None  # noqa: F821
    _resource_path: ClassVar[str]  # Default resource path, subclasses override

    model_config = ConfigDict(
        # populate_by_name=True,
        populate_by_alias=True,
        extra="allow",
        arbitrary_types_allowed=True,
    )

    field_uri: str = Field(..., alias="$uri", exclude=True)

    def __init__(self, **data):
        super().__init__(**data)
        # Initialize resolved cache
        self._resolved_cache = {}

    @computed_field
    def id(self) -> str:
        return self.field_uri.split("/")[-1]

    def __repr__(self):
        return f"<{self.__class__.__name__} {self.id}>"

    def _resolve_and_cache(self, api_ref, cache_key):
        """Resolve an ApiRef and cache the result."""
        if cache_key in self._resolved_cache:
            return self._resolved_cache[cache_key]

        # Check for simple circular reference (a.b.a pattern)
        if hasattr(self, "_immediate_parent") and self._immediate_parent.id == api_ref.id:
            return self._immediate_parent

        # Make API call to resolve, add the parent, cache it
        resolved = api_ref._resolve()
        resolved._immediate_parent = self
        self._resolved_cache[cache_key] = resolved
        return resolved

    def __getattribute__(self, name):
        """Automatically resolve ApiRef objects when accessed."""
        value = super().__getattribute__(name)

        # Only process field attributes - use object.__getattribute__ to avoid recursion
        try:
            model_fields = object.__getattribute__(self, "model_fields")
        except AttributeError:
            return value

        if name not in model_fields:
            return value

        # Handle ApiBaseModel objects (`?expand` data) - set up parent context
        if isinstance(value, ApiBaseModel):
            value._immediate_parent = self
            return value

        # Handle ApiRef objects
        if isinstance(value, ApiRef):
            return self._resolve_and_cache(value, f"{name}:{value.ref}")

        # Handle lists of ApiRef objects
        if isinstance(value, list) and all(isinstance(item, ApiRef) for item in value):
            return [self._resolve_and_cache(item, f"{name}:{item.ref}") for item in value]

        return value

    @classmethod
    def all(cls) -> List["ApiBaseModel"]:
        raise NotImplementedError("Subclasses must implement this method")

    @classmethod
    def get(cls, id: str) -> "ApiBaseModel":
        resp = cls._client.get(f"{cls._api._base_url}{cls._resource_path}/{id}?expand=all")
        return cls.model_validate(resp.json())

    @classmethod
    def where(cls, **kwargs) -> List["ApiBaseModel"]:
        raise NotImplementedError("Subclasses must implement this method")

    @classmethod
    def create(cls, **kwargs):
        raise NotImplementedError("Subclasses must implement this method")

    def update(self, id: str, **kwargs):
        raise NotImplementedError("Subclasses must implement this method")

    def delete(self, id: str):
        raise NotImplementedError("Subclasses must implement this method")

    def save(self):
        raise NotImplementedError("Subclasses must implement this method")
