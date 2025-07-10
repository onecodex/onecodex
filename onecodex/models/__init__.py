"""
One Codex models module.

This module exports the new Pydantic-based models that replace
the Potion-Client based implementation.
"""

# Import base classes
from onecodex.models.base import OneCodexModel, ResourceList

# For backward compatibility, provide OneCodexBase as alias
OneCodexBase = OneCodexModel

# Import all models
from onecodex.models.analyses import (  # noqa
    Analyses,
    Classifications,
    Alignments,
    Panels,
    FunctionalProfiles,
)
from onecodex.models.collection import SampleCollection  # noqa
from onecodex.models.misc import Assets, Jobs, Projects, Tags, Users, Documents  # noqa
from onecodex.models.samples import Samples, Metadata  # noqa

# Import experimental models
from onecodex.models.experimental import (  # noqa
    AnnotationSets,
    Assemblies,
    Genomes,
    Taxa,
)

# Model lookup for API binding
_model_lookup = {}

# Export all models
__all__ = [
    "OneCodexBase",
    "ResourceList",
    "Alignments",
    "Assets",
    "Classifications",
    "Documents",
    "FunctionalProfiles",
    "Jobs",
    "Metadata",
    "Panels",
    "Projects",
    "SampleCollection",
    "Samples",
    "Tags",
    "Users",
    "AnnotationSets",
    "Assemblies",
    "Genomes",
    "Taxa",
]

# Build model lookup (deferred to avoid circular import issues)
import inspect
import sys


def build_model_lookup():
    """Build the model lookup table."""
    global _model_lookup

    def is_oc_class(cls):
        return inspect.isclass(cls) and issubclass(cls, OneCodexModel)

    for name, obj in inspect.getmembers(sys.modules[__name__], is_oc_class):
        if hasattr(obj, "_resource_path"):
            # Get the actual string value, not the descriptor
            try:
                resource_path = obj._resource_path
                if isinstance(resource_path, str) and resource_path:
                    _model_lookup[resource_path] = obj
            except Exception:
                pass  # Skip models without proper _resource_path


# Utility function for error formatting
def pretty_print_error(err_json):
    """Pretty print error messages for the user."""
    # Special case validation errors
    if len(err_json) == 1 and "validationOf" in err_json[0]:
        required_fields = ", ".join(err_json[0]["validationOf"]["required"])
        return f"Validation error. Requires properties: {required_fields}."

    # General error handling
    msg = "; ".join(err.get("message", "") for err in err_json)

    # Fallback
    if not msg:
        msg = "Bad request."
    return msg
