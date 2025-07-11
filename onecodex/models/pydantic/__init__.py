from onecodex.models.pydantic.analysis import (
    Classifications,
    Analyses,
    Alignments,
    FunctionalProfiles,
    Panels,
)
from onecodex.models.pydantic.sample import Metadata, Samples
from onecodex.models.pydantic.misc import Users, Projects, Tags, Jobs, Documents

# Rebuild all models to resolve forward references
Metadata.model_rebuild()
Samples.model_rebuild()
Classifications.model_rebuild()
Analyses.model_rebuild()
Alignments.model_rebuild()
FunctionalProfiles.model_rebuild()
Panels.model_rebuild()
Users.model_rebuild()
Projects.model_rebuild()
Tags.model_rebuild()
Jobs.model_rebuild()
Documents.model_rebuild()

__all__ = [
    "Samples",
    "Metadata",
    "Users",
    "Projects",
    "Tags",
    "Classifications",
    "Alignments",
    "FunctionalProfiles",
    "Panels",
    "Analyses",
    "Jobs",
    "Documents",
]
