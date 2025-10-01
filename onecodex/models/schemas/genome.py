from typing import List, Union

from onecodex.models.base import ApiRef
from onecodex.models.schemas.base import URIModel
from onecodex.models.schemas.misc import JobSchema, TagSchema, UserSchema
from onecodex.models.schemas.sample import ApiV1Visibility, SampleSchema
from onecodex.models.schemas.types import RFC3339Datetime


class AssemblySchema(URIModel):
    created_at: RFC3339Datetime
    filename: str | None
    genome: ApiRef | None
    input_samples: Union[List[SampleSchema], List[ApiRef]] | None
    job: Union[JobSchema, ApiRef] | None
    owner: Union[UserSchema, ApiRef]
    primary_annotation_set: ApiRef | None
    size: int | None
    visibility: ApiV1Visibility


class AnnotationSetSchema(URIModel):
    created_at: RFC3339Datetime
    assembly: Union[AssemblySchema, ApiRef]
    job: Union[JobSchema, ApiRef] | None


class TaxonSchema(URIModel):
    created_at: RFC3339Datetime
    name: str | None
    taxon_id: str
    rank: str | None
    parent: ApiRef | None


class GenomeSchema(URIModel):
    created_at: RFC3339Datetime
    assemblies: Union[List[AssemblySchema], List[ApiRef]]
    description: str | None
    name: str | None
    primary_assembly: Union[AssemblySchema, ApiRef] | None
    tags: Union[List[TagSchema], List[ApiRef]]
    taxon: Union[TaxonSchema, ApiRef]
