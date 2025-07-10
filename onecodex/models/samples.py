"""
Sample model classes extending the generated models.
"""

from typing import Any, Dict, List, Optional

from onecodex.models.base import OneCodexModel, ResourceList
from onecodex.models.mixins import DownloadMixin, UploadMixin, truncate_string
from onecodex.models.generated import SampleSchema, UserSchema, ProjectSchema


class Samples(OneCodexModel, DownloadMixin, UploadMixin):
    """Sample model with One Codex-specific functionality."""
    
    _resource_path = "/api/v1/samples"
    
    # Core fields from SampleSchema
    filename: Optional[str] = None
    size: Optional[int] = None
    created_at: Optional[str] = None
    updated_at: Optional[str] = None
    visibility: Optional[str] = None
    
    # Relationships (will be resolved to model instances)
    # Note: Using Any for now to avoid circular import issues
    project: Optional[Any] = None
    user: Optional[Any] = None
    
    def __repr__(self) -> str:
        """String representation of the sample."""
        filename = getattr(self, 'filename', None) or "(N/A)"
        return f'<{self.__class__.__name__} {self.id}: "{truncate_string(filename, 24)}">'
    
    @classmethod
    def where(cls, *args, **kwargs):
        """Query and retrieve a set of samples.
        
        Args:
            limit: Maximum number of samples to retrieve
            organization: Search all samples within your organization
            public: Search all public samples (limited to 1000 results)
            tags: List of Tags to filter by
            project: Project to filter by
            **kwargs: Additional filter parameters
            
        Returns:
            SampleCollection with samples matching the query
        """
        from onecodex.models.new_collection import SampleCollection
        
        # Handle special parameters
        public = kwargs.pop("public", False)
        organization = kwargs.pop("organization", False)
        tags = kwargs.pop("tags", [])
        project = kwargs.pop("project", None)
        
        # Convert tags and project to appropriate filters
        if tags:
            if not isinstance(tags, list):
                tags = [tags]
            # TODO: Implement tag filtering
            # For now, add as a regular filter
            kwargs["tags"] = [getattr(tag, 'id', tag) for tag in tags]
        
        if project:
            # Convert project to ID if it's a model instance
            project_id = getattr(project, 'id', project)
            kwargs["project"] = project_id
        
        # Set default limit for public searches
        if public and "limit" not in kwargs:
            kwargs["limit"] = 1000
        
        # TODO: Handle organization and public search routes
        # For now, use standard where method
        models = super().where(*args, **kwargs)
        return SampleCollection([m for m in models], cls)


class Metadata(OneCodexModel):
    """Metadata model for samples."""
    
    _resource_path = "/api/v1/metadata"
    
    # Metadata fields
    sample: Optional[Any] = None
    metadata: Optional[Dict[str, Any]] = None