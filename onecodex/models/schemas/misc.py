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
    owner: Optional[Union["UserSchema", ApiRef]]
    permissions: list[str]
    project_name: Optional[str] = Field(
        pattern=PROJECT_HANDLE_REGEX.pattern,
        min_length=PROJECT_HANDLE_MIN_LEN,
        max_length=PROJECT_HANDLE_MAX_LEN,
    )
    public: bool = False


class RepositorySchema(BaseModel):
    url: str = Field(description="The URL of the git repository (https only).")
    tag: Optional[str] = Field(default=None, description="The git tag to use for this repository.")


class JobDependencyRef(BaseModel):
    job: ApiRef
    output_dir: str


class _JobMutableFields(BaseModel):
    """Fields that may be set when creating or updating a Job.

    Optional here so the type is shared between Create and Update schemas;
    Create tightens the required fields.
    """

    name: Optional[str] = None
    script: Optional[str] = None
    image_uri: Optional[str] = None
    description: Optional[str] = None
    cpu: Optional[float] = None
    ram_gb: Optional[float] = None
    storage_gb: Optional[float] = None
    inject_bearer_token: Optional[bool] = None
    repository: Optional[RepositorySchema] = None
    assets: Optional[list[ApiRef]] = None
    dependencies: Optional[list[JobDependencyRef]] = None
    arguments_schema: Optional[list[dict[str, Any]]] = None


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
    job_type: Optional[str] = None


class JobDetails(BaseModel):
    """These are the editable inputs on a Job."""

    model_config = ConfigDict(extra="ignore")

    description: str
    job_type: str
    script: str
    image_uri: str
    cpu: float
    ram_gb: float
    storage_gb: float
    repository: RepositorySchema
    assets: list[ApiRef]
    dependencies: list[JobDependencyRef]
    arguments_schema: list[dict[str, Any]]
    inject_bearer_token: bool
    autorun_on_org_sample_upload: bool


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


class CostSchema(BaseModel):
    amount: int
    currency: str


class AssetUpdateSchema(BaseModel):
    name: str


class JobCreateSchema(_JobMutableFields):
    model_config = ConfigDict(extra="forbid")

    name: str
    script: str
    image_uri: str
    job_type: Optional[str] = None


class JobUpdateSchema(_JobMutableFields):
    model_config = ConfigDict(extra="forbid")

    autorun_on_org_sample_upload: Optional[bool] = None


class AssetSchema(URIModel):
    created_at: RFC3339Datetime
    name: str
    filename: str
    size: Optional[int] = None
    status: AssetStatus
    uploader: Union[UserSchema, ApiRef]
