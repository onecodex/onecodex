from pydantic import BaseModel, ConfigDict, Field


class URIModel(BaseModel):
    model_config = ConfigDict(
        populate_by_alias=True,
        extra="ignore",
        arbitrary_types_allowed=False,
    )
    field_uri: str = Field(alias="$uri")
