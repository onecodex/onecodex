import warnings
from typing import Any, List, Optional, Union, TYPE_CHECKING

from pydantic import Field
from dataclasses import dataclass

from onecodex.exceptions import MethodNotSupported, OneCodexException
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

if TYPE_CHECKING:
    from onecodex.models.sample import Samples
    from onecodex.models.analysis import Analyses


@dataclass(frozen=True)
class DependencyOverride:
    analysis_id: str
    download_path: str | None


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

    def run(
        self,
        sample: "Samples",
        job_args: dict[str, Any] | None = None,
        dependency_overrides: list[DependencyOverride] | None = None,
        populate_default_arguments: bool = True,
    ) -> "Analyses":
        from onecodex.models.analysis import Analyses

        url = f"{self._api._base_url}{self._resource_path}/{self.id}/run"
        payload: dict[str, str | dict | list | bool] = {
            "sample": sample.id,
            "job_args": job_args,
            "populate_default_arguments": populate_default_arguments,
        }
        if dependency_overrides:
            payload["dependencies"] = [
                {"analysis": dep.analysis_id, "download_path": dep.download_path}
                for dep in dependency_overrides
            ]

        resp = self._client.post(url, json=payload)
        resp.raise_for_status()
        if "$ref" not in resp.json():
            raise OneCodexException(f"Invalid response when running job {self.id}")

        analysis_id = resp.json()["$ref"].split("/")[-1]
        return Analyses.get(analysis_id)

    def __repr__(self):
        return f"<{self.__class__.__name__} {self.name} - {self.analysis_type} ({self.id})>"


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
