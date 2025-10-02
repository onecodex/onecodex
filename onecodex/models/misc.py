import warnings
from typing import List, Optional, Union

from pydantic import Field

from onecodex.exceptions import MethodNotSupported
from onecodex.lib.upload import upload_document, upload_asset
from onecodex.models.helpers import truncate_string, ResourceDownloadMixin

from onecodex.models.base import OneCodexBase, ApiRef
from onecodex.models.schemas.misc import (
    TagSchema,
    UserSchema,
    ProjectSchema,
    JobSchema,
    DocumentSchema,
    AssetSchema,
    AssetUpdateSchema,
)


class Tags(OneCodexBase, TagSchema):
    _resource_path = "/api/v1/tags"
    field_uri: Optional[str] = Field(None, alias="$uri")

    def __repr__(self):
        return '<{} {}: "{}">'.format(
            self.__class__.__name__, self.id, truncate_string(self.name, 24)
        )

    def __hash__(self):
        return hash(self.name)

    def save(self):
        # Tags cannot be saved directly this way; instruct user to save via a sample
        raise MethodNotSupported(
            "Tags cannot be saved directly. Instead, append a newly created tag to a sample and save the sample."
        )


class Users(OneCodexBase, UserSchema):
    _resource_path = "/api/v1/users"


class Projects(OneCodexBase, ProjectSchema):
    _resource_path = "/api/v1/projects"
    _allowed_methods = {
        "delete": None,
        "instances_public": None,
    }

    owner: Union[Users, ApiRef]

    @classmethod
    def search_public(cls, *filters, **keyword_filters):
        warnings.warn("Now supported via `.where(..., public=True)`", DeprecationWarning)
        keyword_filters["public"] = True
        keyword_filters["limit"] = 1000
        return cls.where(*filters, **keyword_filters)


class Jobs(OneCodexBase, JobSchema):
    _resource_path = "/api/v1/jobs"


class Documents(OneCodexBase, DocumentSchema, ResourceDownloadMixin):
    _resource_path = "/api/v1/documents"
    _allowed_methods = {
        "delete": None,
    }
    uploader: Union[Users, ApiRef]
    downloaders: List[Union[Users, ApiRef]] = []
    size: Optional[int] = None

    @classmethod
    def upload(cls, file_path, progressbar=None):
        """Upload a series of files to the One Codex server.

        Parameters
        ----------
        file_path : `string`
            A path to a file on the system.
        progressbar : `click.progressbar`, optional
            If passed, display a progress bar using Click.

        Returns
        -------
        A `Documents` object upon successful upload. None if the upload failed.
        """
        doc_id = upload_document(file_path, cls, progressbar=progressbar)
        return cls.get(doc_id)


# NOTE: this model is only accessible via `X-OneCodex-Api-Experimental` as of 10/2/2025
class Assets(OneCodexBase, AssetSchema, ResourceDownloadMixin):
    _resource_path = "/api/v1/assets"
    _allowed_methods = {
        "delete": None,
        "update": AssetUpdateSchema,
    }
    uploaded_by: Union[Users, ApiRef]

    @classmethod
    def upload(cls, file_path, progressbar=None, name=None):
        """Upload a file to an asset.

        Parameters
        ----------
        file_path : `string`
            A path to a file on the system.
        progressbar : `click.progressbar`, optional
            If passed, display a progress bar using Click.
        name : `string`, optional
            If passed, name is sent with upload request and is associated with asset.

        Returns
        -------
        An `Assets` object upon successful upload. None if the upload failed.
        """
        asset_id = upload_asset(file_path, cls, progressbar=progressbar, name=name)

        return cls.get(asset_id)
