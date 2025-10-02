from typing import List, Optional, Union

from onecodex.models.base import ApiRef
from onecodex.models.schemas.base import URIModel
from onecodex.models.schemas.misc import JobSchema, TagSchema, UserSchema
from onecodex.models.schemas.sample import ApiV1Visibility, SampleSchema
from onecodex.models.schemas.types import RFC3339Datetime

# NOTE: these models are only accessible via `X-OneCodex-Api-Experimental` as of 10/2/2025


class AssemblySchema(URIModel):
    created_at: RFC3339Datetime
    filename: Optional[str]
    genome: Optional[ApiRef]
    input_samples: Optional[Union[List[SampleSchema], List[ApiRef]]]
    job: Optional[Union[JobSchema, ApiRef]]
    owner: Union[UserSchema, ApiRef]
    primary_annotation_set: Optional[ApiRef]
    size: Optional[int]
    visibility: ApiV1Visibility


class AnnotationSetSchema(URIModel):
    created_at: RFC3339Datetime
    assembly: Union[AssemblySchema, ApiRef]
    job: Optional[Union[JobSchema, ApiRef]]


class TaxonSchema(URIModel):
    created_at: RFC3339Datetime
    name: Optional[str]
    taxon_id: str
    rank: Optional[str]
    parent: Optional[ApiRef]


class GenomeSchema(URIModel):
    created_at: RFC3339Datetime
    assemblies: Union[List[AssemblySchema], List[ApiRef]]
    description: Optional[str]
    name: Optional[str]
    primary_assembly: Optional[Union[AssemblySchema, ApiRef]]
    tags: Union[List[TagSchema], List[ApiRef]]
    taxon: Union[TaxonSchema, ApiRef]
