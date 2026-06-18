import warnings
from datetime import datetime
from typing import Any, List, Optional, Union, TYPE_CHECKING
from typing_extensions import Self

from onecodex.models.base import UNSET
from onecodex.models.filters import (
    DatetimeFilter,
    EqStrFilter,
    NumFilter,
    RefFilter,
    StrFilter,
)

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
    JobCreateSchema,
    JobUpdateSchema,
    JobDetails,
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

    @classmethod
    def where(
        cls,
        *filters: str | dict,
        sort: str | list[str] | None = None,
        limit: int | None = None,
        public: bool = False,
        filter: Any = None,
        name: str | StrFilter = UNSET,
    ) -> list[Self]:
        """Query tags.

        Tags are short labels attached to samples for organizing and
        filtering them.

        Examples
        --------
        Find a tag by exact name::

            tag = ocx.Tags.where(name="trimmed")[0]

        Find tags containing a substring::

            ocx.Tags.where(name={"$icontains": "qc"})

        See :meth:`OneCodexBase.where` for the full operator reference.
        """
        return super().where(
            *filters,
            sort=sort,
            limit=limit,
            public=public,
            filter=filter,
            name=name,
        )

    def save(self):
        # Tags cannot be saved directly this way; instruct user to save via a sample
        raise MethodNotSupported(
            "Tags cannot be saved directly. Instead, append a newly created tag to a sample and save the sample."
        )


class Users(OneCodexBase, UserSchema):
    _resource_path = "/api/v1/users"

    @classmethod
    def where(
        cls,
        *filters: str | dict,
        sort: str | list[str] | None = None,
        limit: int | None = None,
        public: bool = False,
        filter: Any = None,
        email: str | EqStrFilter | None = UNSET,
    ) -> list[Self]:
        """Query users.

        Users are account holders in your organization. The only filterable
        field is ``email``, and the API only supports exact-match (``$eq``)
        — substring operators are not accepted server-side.

        Example
        -------
        Find a user by exact email::

            user = ocx.Users.where(email="alice@example.com")[0]

        See :meth:`OneCodexBase.where` for the full operator reference.
        """
        return super().where(
            *filters,
            sort=sort,
            limit=limit,
            public=public,
            filter=filter,
            email=email,
        )


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
        return cls.where(*filters, public=True, limit=1000, **keyword_filters)

    @classmethod
    def where(
        cls,
        *filters: str | dict,
        sort: str | list[str] | None = None,
        limit: int | None = None,
        public: bool = False,
        filter: Any = None,
        name: str | StrFilter | None = UNSET,
        project_name: str | StrFilter | None = UNSET,
        description: str | StrFilter | None = UNSET,
        external_id: str | StrFilter | None = UNSET,
        owner: Users | str | RefFilter | None = UNSET,
    ) -> list[Self]:
        """Query projects.

        Projects group related samples. Set ``public=True`` to search across
        all public projects on One Codex.

        Examples
        --------
        Find a project by exact name::

            proj = ocx.Projects.where(name="HMP Phase II")[0]

        Find public projects containing a keyword::

            ocx.Projects.where(name={"$icontains": "gut"}, public=True)

        See :meth:`OneCodexBase.where` for the full operator reference.
        """
        return super().where(
            *filters,
            sort=sort,
            limit=limit,
            public=public,
            filter=filter,
            name=name,
            project_name=project_name,
            description=description,
            external_id=external_id,
            owner=owner,
        )


class Jobs(OneCodexBase, JobSchema):
    _resource_path = "/api/v1/jobs"
    _allowed_methods = {
        "create": JobCreateSchema,
        "update": JobUpdateSchema,
    }

    @classmethod
    def where(
        cls,
        *filters: str | dict,
        sort: str | list[str] | None = None,
        limit: int | None = None,
        public: bool = False,
        filter: Any = None,
        created_at: datetime | DatetimeFilter = UNSET,
        name: str | StrFilter = UNSET,
        analysis_type: str | StrFilter = UNSET,
    ) -> list[Self]:
        """Query jobs.

        Jobs are runnable workflows — both built-in (e.g. ``classification``)
        and custom (Nextflow pipelines, shell scripts). Use ``public=True``
        to include One Codex's built-in jobs and any other public custom
        jobs in addition to your own.

        Examples
        --------
        Find your custom jobs by name::

            ocx.Jobs.where(name={"$icontains": "amplicon"})

        Find all classification-type jobs (yours + public)::

            ocx.Jobs.where(analysis_type="classification", public=True)

        See :meth:`OneCodexBase.where` for the full operator reference.
        """
        return super().where(
            *filters,
            sort=sort,
            limit=limit,
            public=public,
            filter=filter,
            created_at=created_at,
            name=name,
            analysis_type=analysis_type,
        )

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
                {"analysis": dep.analysis.id, "download_path": dep.download_path}
                for dep in dependency_overrides
            ]

        resp = self._client.post(url, json=payload)
        if not resp.ok:
            try:
                body = resp.json()
            except ValueError:
                body = None
            detail = None
            if isinstance(body, dict):
                detail = body.get("message") or body.get("msg")
            if not detail:
                detail = (resp.text or "").strip() or None
            msg = f"Job run failed ({resp.status_code})"
            if detail:
                msg = f"{msg}: {detail}"
            raise OneCodexException(msg)
        if "$ref" not in resp.json():
            raise OneCodexException(f"Invalid response when running job {self.id}")

        analysis_id = resp.json()["$ref"].split("/")[-1]
        return Analyses.get(analysis_id)

    def details(self) -> JobDetails:
        """Fetch the job's detail fields and return them as a `JobDetails`.

        Includes script, image_uri, cpu, ram_gb, storage_gb, repository, assets,
        dependencies, and arguments_schema. Only available for user-created jobs the
        current user has permission to view the details of.
        """
        url = f"{self._api._base_url}{self._resource_path}/{self.id}/details"
        resp = self._client.get(url)
        resp.raise_for_status()
        return JobDetails.model_validate(resp.json())

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
    def where(
        cls,
        *filters: str | dict,
        sort: str | list[str] | None = None,
        limit: int | None = None,
        public: bool = False,
        filter: Any = None,
        created_at: datetime | DatetimeFilter = UNSET,
        filename: str | StrFilter = UNSET,
        size: int | NumFilter | None = UNSET,
        uploader: Users | str | RefFilter = UNSET,
    ) -> list[Self]:
        """Query documents.

        Documents are arbitrary files (PDFs, spreadsheets, etc.) uploaded to
        your account, optionally shared with other users.

        Examples
        --------
        Find documents uploaded since a date::

            from datetime import datetime, timedelta, timezone
            since = (datetime.now(timezone.utc) - timedelta(days=7)).isoformat()
            ocx.Documents.where(created_at={"$gte": since})

        Find documents larger than 10 MB::

            ocx.Documents.where(size={"$gte": 10 * 1024 * 1024})

        Find documents uploaded by a specific user::

            ocx.Documents.where(uploader=user)

        See :meth:`OneCodexBase.where` for the full operator reference.
        """
        return super().where(
            *filters,
            sort=sort,
            limit=limit,
            public=public,
            filter=filter,
            created_at=created_at,
            filename=filename,
            size=size,
            uploader=uploader,
        )

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
    def where(
        cls,
        *filters: str | dict,
        sort: str | list[str] | None = None,
        limit: int | None = None,
        public: bool = False,
        filter: Any = None,
        created_at: datetime | DatetimeFilter = UNSET,
        name: str | StrFilter = UNSET,
        filename: str | StrFilter = UNSET,
        size: int | NumFilter | None = UNSET,
        status: str | StrFilter = UNSET,
        uploader: Users | str | RefFilter = UNSET,
    ) -> list[Self]:
        """Query assets.

        Assets are reusable inputs to custom jobs — reference databases,
        index files, anything you want available at job run time. Shared
        with everyone in your organization.

        Examples
        --------
        Find assets by name::

            ocx.Assets.where(name={"$icontains": "kraken"})

        Find only assets ready to use (``status="available"``)::

            ocx.Assets.where(status="available")

        See :meth:`OneCodexBase.where` for the full operator reference.
        """
        return super().where(
            *filters,
            sort=sort,
            limit=limit,
            public=public,
            filter=filter,
            created_at=created_at,
            name=name,
            filename=filename,
            size=size,
            status=status,
            uploader=uploader,
        )

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
