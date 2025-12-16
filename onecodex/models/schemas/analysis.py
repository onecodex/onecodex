from __future__ import annotations
from typing import Any, Optional, Union
from pydantic import Field, ConfigDict

from onecodex.models.base import ApiRef
from onecodex.models.schemas.base import URIModel
from onecodex.models.schemas.types import RFC3339Datetime

from onecodex.models.schemas.misc import CostSchema, JobSchema
from onecodex.models.schemas.sample import SampleSchema


class BaseAnalysisSchema(URIModel):
    model_config = ConfigDict(extra="ignore")

    created_at: RFC3339Datetime
    complete: bool = False
    error_msg: Optional[str] = Field(
        default=None,
        description="The error message from the analysis, if it failed.",
        max_length=255,
    )
    job: Union["JobSchema", ApiRef] = Field(  # noqa: F821
        description='A reference to the versioned job underlying the analysis, e.g., `{"$ref": "/api/v1/jobs/d512cb556241440f"}`.',
    )
    job_args: dict[str, Any] = Field(
        default={},
        description="The arguments passed into this analysis (can be `null`).",
    )
    sample: Union["SampleSchema", ApiRef] = Field(  # noqa: F821
        description='A reference to the sample underlying the analysis, e.g., `{"$ref": "/api/v1/sample/0ee172af60e84f61"}`.',
    )
    success: Optional[bool] = False

    cost: Optional[CostSchema]
    draft: bool
    dependencies: Union[list[BaseAnalysisSchema], list[ApiRef]]


class AnalysisSchema(BaseAnalysisSchema):
    analysis_type: str


class AlignmentSchema(BaseAnalysisSchema):
    pass


class ClassificationSchema(BaseAnalysisSchema):
    pass


class FunctionalRunSchema(BaseAnalysisSchema):
    pass


class PanelSchema(BaseAnalysisSchema):
    pass
