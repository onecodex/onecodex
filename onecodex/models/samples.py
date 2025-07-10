"""
Sample model classes extending the generated models.
"""

from typing import Any, Dict, List, Optional

from onecodex.models.base import OneCodexModel, ResourceList
from onecodex.models.generated import SampleSchema, UserSchema, ProjectSchema


class Samples(OneCodexModel):
    """Sample model with One Codex-specific functionality."""
    
    _resource_path = "/api/v1/samples"
    
    # Core fields from SampleSchema
    filename: Optional[str] = None
    size: Optional[int] = None
    created_at: Optional[str] = None
    updated_at: Optional[str] = None
    visibility: Optional[str] = None
    
    # Relationships (will be resolved to model instances)
    project: Optional["Projects"] = None
    user: Optional["Users"] = None
    
    @classmethod
    def where(cls, *args, **kwargs):
        """Override where to return SampleCollection."""
        from onecodex.models.collection import SampleCollection
        
        models = super().where(*args, **kwargs)
        return SampleCollection([m for m in models], cls)


class Metadata(OneCodexModel):
    """Metadata model for samples."""
    
    _resource_path = "/api/v1/metadata"
    
    # Metadata fields
    sample: Optional[Samples] = None
    metadata: Optional[Dict[str, Any]] = None


# Forward reference resolution
Samples.model_rebuild()