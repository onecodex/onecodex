import re

from typing import Any, List, Optional, Union
from pydantic import Field, ConfigDict, BaseModel

from onecodex.models.base import ApiRef
from onecodex.models.schemas.base import URIModel
from onecodex.models.schemas.types import RFC3339Datetime


# FIXME: Make these literals or enums
AssetStatus = str


class TagSchema(URIModel):
    name: str = Field(
        description="The tag label or name. Name must be 30 characters or fewer.",
        max_length=30,
    )


class UserSchema(URIModel):
    email: Optional[str]


PROJECT_HANDLE_MIN_LEN = 3
PROJECT_HANDLE_MAX_LEN = 15
PROJECT_HANDLE_REGEX = re.compile("^[a-zA-Z0-9_-]{3,15}$")


class ProjectSchema(URIModel):
    description: Optional[str] = None
    external_id: Optional[str] = Field(None, json_schema_extra={"example": "a1b2c3d4e5f67890"})
    name: Optional[str] = None
    owner: Union["UserSchema", ApiRef]
    permissions: list[str]
    project_name: Optional[str] = Field(
        pattern=PROJECT_HANDLE_REGEX.pattern,
        min_length=PROJECT_HANDLE_MIN_LEN,
        max_length=PROJECT_HANDLE_MAX_LEN,
    )
    public: bool = False


class JobSchema(URIModel):
    created_at: RFC3339Datetime
    name: str = Field(
        description="The name of the job (this is displayed in the dropdown on the analysis page of the One Codex web application)."
    )
    job_args_schema: dict[str, Any] = Field(
        default_factory=dict,
        description="The JSON schema for the arguments taken by the job (can be an empty object, i.e., `{}`).",
    )
    analysis_type: str = Field(
        description="The type of analysis. See [the analysis resource](/api/analysis-resource) for more details."
    )
    public: bool = Field(
        description="Whether the job is publicly available. For most jobs this will be `true`. Custom, private jobs are also available, and will only be visible to users whose samples (or samples shared with them) have been analyzed using that job."
    )


class DocumentSchema(URIModel):
    # Ignore `description` and `valid_until` fields
    model_config = ConfigDict(extra="ignore")

    created_at: RFC3339Datetime
    filename: str = Field(description='The document filename (e.g., "report.pdf")')
    size: Optional[int] = Field(description="The size of the document in bytes.", ge=1)
    uploader: Union["UserSchema", ApiRef] = Field(
        description='A reference to the user that uploaded and owns the document, e.g., `{"$ref": "/api/v1/users/5891ee65711c4d5e"}`. Only owners can modify a document (with some exceptions in the case of organization accounts configured for multiple users - please [contact us](mailto:support@onecodex.com) if you\'d like to discuss this use case).',
    )
    downloaders: Union[List["UserSchema"], List[ApiRef]] = Field(
        description="An (optionally empty) array of references to users that the document has been shared with. These users are able to download (but not modify) the document.",
    )


class FileDetailSchema(BaseModel):
    filename: str
    filepath: str
    size: int
    url: str


# NOTE: these Asset models are only accessible via `X-OneCodex-Api-Experimental` as of 10/2/2025


class CostSchema(BaseModel):
    amount: float
    currency: str


class AssetUpdateSchema(BaseModel):
    name: str


# NOTE: this model is unstable/experimental and likely to change
class AssetSchema(URIModel):
    created_at: RFC3339Datetime
    name: str
    filename: str
    s3_uri: Optional[str]
    status: AssetStatus
    organization_id: int
    uploaded_by: Union[UserSchema, ApiRef]
    uuid: str
