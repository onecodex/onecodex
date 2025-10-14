from onecodex.models.analysis import (
    Classifications,
    Analyses,
    Alignments,
    FunctionalProfiles,
    Panels,
)
from onecodex.models.sample import Metadata, Samples
from onecodex.models.collection import SampleCollection
from onecodex.models.misc import Users, Projects, Tags, Jobs, Documents, Assets
from onecodex.models.genome import AnnotationSets, Assemblies, Genomes, Taxa

# Model registry for dynamic resolution based on resource paths
_MODEL_REGISTRY = {}


def register_model(cls):
    """Register a model class in the global registry."""
    if hasattr(cls, "_resource_path"):
        _MODEL_REGISTRY[cls._resource_path] = cls
    return cls


def get_model_class(uri):
    """Get the model class for a given URI."""
    for resource_path, model_class in _MODEL_REGISTRY.items():
        if resource_path in uri:
            return model_class
    return None


# Register all model classes
register_model(Samples)
register_model(Metadata)
register_model(Classifications)
register_model(Analyses)
register_model(Alignments)
register_model(FunctionalProfiles)
register_model(Panels)
register_model(Users)
register_model(Projects)
register_model(Tags)
register_model(Jobs)
register_model(Documents)
register_model(Assets)
register_model(AnnotationSets)
register_model(Assemblies)
register_model(Genomes)
register_model(Taxa)

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
Assets.model_rebuild()
AnnotationSets.model_rebuild()
Assemblies.model_rebuild()
Genomes.model_rebuild()
Taxa.model_rebuild()

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
    "Assets",
    "AnnotationSets",
    "Assemblies",
    "Genomes",
    "Taxa",
    "SampleCollection",
]
