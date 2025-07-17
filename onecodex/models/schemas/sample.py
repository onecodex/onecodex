from typing import Any, Optional, List, Union
from pydantic import Field, ConfigDict, BaseModel

from onecodex.models.base import ApiRef
from onecodex.models.schemas.base import URIModel
from onecodex.models.schemas.types import RFC3339Datetime

# FIXME: Make these literals or enums
SequencingPlatform = str
SequencingSampleType = str
SequencingLibraryType = str
UserFacingDataStatus = str
ApiV1Visibility = str


class SampleSchema(URIModel):
    created_at: RFC3339Datetime
    error_msg: Optional[str] = Field(
        None,
        description="An error message if the sample failed to upload, import, or validate.",
    )
    filename: Optional[str] = Field(
        None,
        description="The sample filename.",
    )
    metadata: Union["MetadataSchema", ApiRef] = Field(
        description="A metadata object.",
    )
    owner: Union["UserSchema", ApiRef] = Field(  # noqa: F821
        description="The owner of the sample.",
    )
    primary_classification: Optional[Union["ClassificationSchema", ApiRef]] = (  # noqa: F821
        Field(
            description="A reference to a Classification for the sample. This will typically be the [One Codex Database](https://docs.onecodex.com/en/articles/3761205-one-codex-database) or [Targeted Loci Database](https://docs.onecodex.com/en/articles/3754219-targeted-loci-database) results as appropriate. Note that samples will not have a `primary_classification` while they are still importing or being uploaded.",
        )
    )
    project: Optional[Union["ProjectSchema", ApiRef]] = Field(  # noqa: F821
        None,
        description="The project the sample belongs to (optional).",
    )
    size: Optional[int] = Field(
        None,
        description="The size of the uploaded file in bytes.",
    )
    tags: List[Union["TagSchema", ApiRef]] = Field(  # noqa: F821
        default_factory=list,
        description="An (optionally empty) array of references to Tags describing the sample. Tags are an additional unstructured organizational tool that complement Projects and Metadata records.",
    )
    status: UserFacingDataStatus = Field(
        description="The status of the sample.",
    )
    visibility: ApiV1Visibility = Field(
        description="The visibility of the sample (affects who can view the sample and its analyses).",
    )
    updated_at: Optional[RFC3339Datetime] = None


class SampleUpdateSchema(BaseModel):
    model_config = ConfigDict(extra="ignore")

    project: Optional[ApiRef] = None  # noqa: F821
    tags: List[ApiRef] = Field(default_factory=list)  # noqa: F821
    visibility: ApiV1Visibility = Field(
        default="private",
        description="The visibility of the sample (affects who can view the sample and its analyses).",
    )


class _MetadataOptionalFields(BaseModel):
    custom: Optional[dict[str, Any]] = Field(
        default={},
        description='Arbitrary metadata is supported as part of a custom object. custom has two constraints: (1) it must have a depth of one (i.e., no nested records); and (2) only strings, numbers, boolean, and null values are supported as values. Example: `{"lab_tech": "Linus Pauling", "amplicon_scheme": "V3-V4"}`',
    )
    date_collected: Optional[RFC3339Datetime] = Field(
        default=None, description="Timestamp for when the sample was collected."
    )
    date_sequenced: Optional[RFC3339Datetime] = Field(
        default=None, description="Timestamp for when the sample was sequenced."
    )
    description: Optional[str] = None
    external_sample_id: Optional[str] = Field(
        default=None,
        description="An arbitrary external sample ID, e.g., an ID in a LIMS. Up to 60 characters.",
        json_schema_extra={"example": "a1b2c3d4e5f67890"},
    )
    library_type: Optional[SequencingLibraryType] = Field(
        default=None, description="An enum with the sample library type."
    )
    location_lat: Optional[float] = Field(
        default=None,
        ge=-90.0,
        le=90.0,
        description="The latitude `(-90.0-90.0)` of the sample location. By convention, we recommend using this for the location in which the physical specimen was collected.",
    )
    location_lon: Optional[float] = Field(
        default=None,
        ge=-180,
        le=180,
        description="The longitude `(-180.0-180.0)` of the sample location.",
    )
    location_string: Optional[str] = Field(default=None, max_length=255)
    name: Optional[str] = Field(default=None, max_length=255)
    platform: Optional[SequencingPlatform] = Field(
        default=None, description="An enum with the name of the sequencing platform."
    )
    sample_type: Optional[SequencingSampleType] = Field(
        default=None, description="An enum with the sample type."
    )


class MetadataSchema(_MetadataOptionalFields, URIModel):
    sample: Union[SampleSchema, ApiRef] = Field(  # noqa: F821 # type: ignore[unresolved-reference]
        description="The sample the metadata belongs to.",
    )
    starred: bool = Field(
        default=False,
        description="Whether the sample has been starred by the user within the One Codex web application.",
    )
    updated_at: Optional[RFC3339Datetime] = None


class MetadataPatchSchema(_MetadataOptionalFields):
    starred: Optional[bool] = Field(
        default=False,
        description="Whether the sample has been starred by the user within the One Codex web application.",
    )
