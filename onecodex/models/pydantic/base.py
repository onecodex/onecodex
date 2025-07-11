from typing import ClassVar, Optional, List

from pydantic import BaseModel as PydanticBaseModel
from pydantic import ConfigDict, computed_field, Field, model_validator


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

    @computed_field
    def id(self) -> str:
        return self.field_uri.split("/")[-1]

    def __repr__(self):
        return f"<{self.__class__.__name__} {self.id}>"

    @model_validator(mode="before")
    @classmethod
    def resolve_references(cls, data):
        """Convert $ref dicts to ApiRef objects."""
        if not isinstance(data, dict):
            return data

        current_uri = data.get("$uri")
        print(f"PROCESSING {cls.__name__}: {current_uri}")

        # Process each field that might contain references
        for field_name, field_value in list(data.items()):
            if field_name in ("$uri", "created_at", "updated_at"):
                continue

            data[field_name] = cls._convert_refs_to_objects(field_value)

        return data

    @classmethod
    def _convert_refs_to_objects(cls, value):
        """Convert $ref dicts to ApiRef objects."""
        if isinstance(value, dict) and len(value) == 1 and "$ref" in value:
            # Convert single $ref dict to ApiRef object
            return ApiRef.model_validate(value)
        elif isinstance(value, list):
            # Convert any $ref dicts in lists to ApiRef objects
            return [
                ApiRef.model_validate(item)
                if isinstance(item, dict) and len(item) == 1 and "$ref" in item
                else item
                for item in value
            ]
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
