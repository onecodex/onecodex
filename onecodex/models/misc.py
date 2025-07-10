import warnings

from onecodex.exceptions import OneCodexException, MethodNotSupported
from onecodex.models.base import OneCodexModel
from onecodex.models.mixins import DownloadMixin, UploadMixin, truncate_string


class Tags(OneCodexModel):
    _resource_path = "/api/v1/tags"
    
    # Core tag fields
    name: str = None
    color: str = None
    created_at: str = None

    def __init__(self, **data):
        # Handle special case where Tags can be looked up by name
        if "name" in data and "id" not in data:
            name = data["name"]
            sample = data.pop("sample", None)
            
            # Try to lookup existing tag
            if hasattr(self.__class__, '_api') and self.__class__._api:
                try:
                    results = self.__class__.where(name=name)
                    if len(results) == 1:
                        # Use existing tag
                        super().__init__(**results[0].model_dump())
                        return
                    elif len(results) > 1:
                        raise OneCodexException("Multiple matches found for given criteria")
                except Exception:
                    pass  # Fall through to create new tag
            
            # Create new tag
            if sample:
                data["sample"] = sample
        
        super().__init__(**data)

    def __repr__(self):
        name = getattr(self, 'name', 'Unknown')
        return f'<{self.__class__.__name__} {self.id}: "{truncate_string(name, 24)}">'

    def __hash__(self):
        return hash(getattr(self, 'name', ''))

    def save(self):
        # Tags cannot be saved directly this way; instruct user to save via a sample
        raise MethodNotSupported(
            "Tags cannot be saved directly. Instead, append a newly created tag to a sample and save the sample."
        )


class Users(OneCodexModel):
    _resource_path = "/api/v1/users"
    
    # Core user fields
    name: str = None
    email: str = None
    created_at: str = None
    updated_at: str = None


class Projects(OneCodexModel):
    _resource_path = "/api/v1/projects"
    
    # Core project fields
    name: str = None
    description: str = None
    created_at: str = None
    updated_at: str = None
    visibility: str = None

    @classmethod
    def search_public(cls, *filters, **keyword_filters):
        warnings.warn("Now supported via `.where(..., public=True)`", DeprecationWarning)
        keyword_filters["public"] = True
        keyword_filters["limit"] = 1000
        return cls.where(*filters, **keyword_filters)


class Jobs(OneCodexModel):
    _resource_path = "/api/v1/jobs"
    
    # Core job fields
    status: str = None
    created_at: str = None
    updated_at: str = None
    error_msg: str = None


class Documents(OneCodexModel, DownloadMixin, UploadMixin):
    _resource_path = "/api/v1/documents"
    
    # Core document fields
    filename: str = None
    size: int = None
    content_type: str = None
    created_at: str = None
    updated_at: str = None

    @classmethod
    def upload(cls, file_path, progressbar=None):
        """Upload a file to the One Codex server.

        Args:
            file_path: Path to a file on the system.
            progressbar: Optional Click progressbar for progress display.

        Returns:
            A Documents object upon successful upload. None if the upload failed.
        """
        # Use the UploadMixin functionality
        return super().upload(file_path, progressbar=progressbar)


class Assets(OneCodexModel, DownloadMixin, UploadMixin):
    _resource_path = "/api/v1_experimental/assets"
    
    # Core asset fields
    filename: str = None
    size: int = None
    content_type: str = None
    created_at: str = None

    @classmethod
    def upload(cls, file_path, progressbar=None, name=None):
        """Upload a file to an asset.

        Args:
            file_path: Path to a file on the system.
            progressbar: Optional Click progressbar for progress display.
            name: Optional name to associate with the asset.

        Returns:
            An Assets object upon successful upload. None if the upload failed.
        """
        # Use the UploadMixin functionality with additional name parameter
        return super().upload(file_path, progressbar=progressbar, name=name)
