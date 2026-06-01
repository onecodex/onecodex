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
    from onecodex.models.analysis import Analyses

    analysis: Analyses
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
        dependency_overrides: "list[Analyses | DependencyOverride] | None" = None,
        populate_default_arguments: bool = True,
    ) -> "Analyses":
        """Run this job against a sample and return the resulting analysis.

        Parameters
        ----------
        sample : `Samples`
            The sample to run this job against.
        job_args : `dict`, optional
            Job arguments, validated server-side against the job's
            ``job_args_schema``. Values are coerced from strings per the schema
            (e.g. an ``integer`` field accepts ``"20"``). Passing arguments that
            do not match the schema raises `OneCodexException` with the
            server's validation message.
        dependency_overrides : `list[Analyses | DependencyOverride]`, optional
            Override one or more of the job's declared dependencies with
            existing analyses. Pass an `Analyses` directly to swap in the
            parent run while keeping the job-declared (or default) output
            path. Wrap in `DependencyOverride` only when you also need to set
            a custom ``download_path``.
        populate_default_arguments : `bool`, default True
            If True, the server fills in defaults for any args not provided in
            ``job_args``. Set False to require every arg to be supplied
            explicitly.

        Returns
        -------
        `Analyses`
            The newly-created analysis. Note: if an analysis with identical
            inputs already exists, the server may return the existing one
            rather than starting a new run.

        Raises
        ------
        OneCodexException
            If the server rejects the request (e.g. invalid ``job_args``,
            permission denied, or unknown sample/dependency). The exception
            message contains the server's explanation.

        Examples
        --------
        Basic usage::

            job = ocx.Jobs.get("47c4fe23588640a9")
            sample = ocx.Samples.get("7428cca4a3a04a8e")
            analysis = job.run(sample)

        With job arguments::

            analysis = job.run(sample, {"min_quality": "20", "adapter": "AGATC"})

        Inspect the job's argument schema first::

            job.job_args_schema  # JSON schema describing accepted args

        Override a dependency to reuse a prior analysis::

            prior = ocx.Analyses.get("abc123def4567890")
            analysis = job.run(sample, dependency_overrides=[prior])

        Override with a custom output path::

            from onecodex.models.misc import DependencyOverride

            analysis = job.run(
                sample,
                dependency_overrides=[
                    DependencyOverride(analysis=prior, download_path="results/parent.tsv"),
                ],
            )
        """
        from onecodex.models.analysis import Analyses

        url = f"{self._api._base_url}{self._resource_path}/{self.id}/run"
        payload: dict[str, Any] = {
            "sample": sample.id,
            "job_args": job_args,
            "populate_default_arguments": populate_default_arguments,
        }
        if dependency_overrides:
            normalized = [
                dep
                if isinstance(dep, DependencyOverride)
                else DependencyOverride(analysis=dep, download_path=None)
                for dep in dependency_overrides
            ]
            payload["dependencies"] = [
                {"analysis": dep.analysis.id, "download_path": dep.download_path}
                for dep in normalized
            ]

        resp = self._client.post(url, json=payload)
        if not resp.ok:
            try:
                message = resp.json().get("message")
            except ValueError:
                message = None
            raise OneCodexException(
                message or f"Error running job {self.id}: HTTP {resp.status_code}"
            )
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


class Assets(OneCodexBase, AssetSchema, ResourceDownloadMixin):
    _resource_path = "/api/v1/assets"
    _allowed_methods = {
        "delete": None,
        "update": AssetUpdateSchema,
    }
    uploader: Union[Users, ApiRef]

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
