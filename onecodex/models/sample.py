from __future__ import annotations

from typing import Any, TYPE_CHECKING, Optional, List, Union
from datetime import datetime

if TYPE_CHECKING:
    from onecodex.models.collection import SampleCollection

from requests.exceptions import HTTPError

from typing_extensions import Self

from onecodex.exceptions import OneCodexException
from onecodex.models.base import UNSET, _drop_unset
from onecodex.models.filters import (
    BoolFilter,
    DatetimeFilter,
    NumFilter,
    RefFilter,
    StrFilter,
)
from onecodex.models.helpers import ResourceDownloadMixin
from onecodex.models.base import OneCodexBase, ApiRef
from onecodex.models.analysis import Classifications

from onecodex.models.helpers import truncate_string
from onecodex.models.schemas import (
    SampleSchema,
    SampleUpdateSchema,
    MetadataSchema,
    MetadataPatchSchema,
)
from onecodex.lib.upload import upload_sequence

from pydantic import field_validator

from onecodex.models.misc import Users, Projects, Tags


def get_project(project):
    """
    Get the actual project instance if the argument is not None and not already a Project.

    Raises an exception if the project can't be found.
    """
    if not isinstance(project, Projects) and project is not None:
        project_search = Projects.get(project)
        if not project_search:
            project_search = Projects.where(name=project)
        if not project_search:
            try:
                project_search = Projects.where(project_name=project)
            except HTTPError:
                project_search = None
        if not project_search:
            raise OneCodexException(
                "Project {} does not exist. Please create the project in One Codex then try again.".format(
                    project
                )
            )

        if isinstance(project_search, list):
            return project_search[0]
        elif isinstance(project_search, Projects):
            return project_search

    return project


class _MetadataSchema(MetadataSchema):
    # Use ApiRef for the circular reference
    sample: Union["Samples", ApiRef]
    updated_at: Optional[datetime] = None
    date_collected: Optional[datetime] = None
    date_sequenced: Optional[datetime] = None

    @field_validator("platform", "library_type", "sample_type")
    @classmethod
    def convert_enum_to_str(cls, v):
        if v is None:
            return v
        return str(v.value) if hasattr(v, "value") else v


class Metadata(OneCodexBase, _MetadataSchema):
    _resource_path = "/api/v1/metadata"
    _allowed_methods = {
        "update": MetadataPatchSchema,
    }
    _use_cursor_pagination = True

    @classmethod
    def where(
        cls,
        *filters: str | dict,
        sort: str | list[str] | None = None,
        limit: int | None = None,
        public: bool = False,
        filter: Any = None,
        updated_at: datetime | DatetimeFilter | None = UNSET,
        date_collected: datetime | DatetimeFilter | None = UNSET,
        date_sequenced: datetime | DatetimeFilter | None = UNSET,
        description: str | StrFilter | None = UNSET,
        external_sample_id: str | StrFilter | None = UNSET,
        library_type: str | StrFilter | None = UNSET,
        location_lat: float | NumFilter | None = UNSET,
        location_lon: float | NumFilter | None = UNSET,
        location_string: str | StrFilter | None = UNSET,
        name: str | StrFilter | None = UNSET,
        platform: str | StrFilter | None = UNSET,
        sample_type: str | StrFilter | None = UNSET,
        starred: bool | BoolFilter = UNSET,
        sample: Samples | str | RefFilter = UNSET,
    ) -> list[Self]:
        """Query sample metadata.

        Metadata holds per-sample structured fields — sequencing platform,
        sample type, collection date, lat/lon, plus an arbitrary ``custom``
        dict. Only metadata for samples you own is queryable; public and
        organization-shared samples cannot be searched by metadata field.

        Examples
        --------
        Find metadata by collection date::

            from datetime import datetime, timedelta, timezone
            cutoff = (datetime.now(timezone.utc) - timedelta(days=30)).isoformat()
            ocx.Metadata.where(date_collected={"$gte": cutoff})

        Find by external (e.g. LIMS) ID::

            md = ocx.Metadata.where(external_sample_id="LIMS-12345")[0]

        ``custom`` is not server-filterable. Fetch a broader set and filter
        in Python::

            ocx.Samples.where(...).filter(lambda s: s.metadata.custom.get("lab") == "X")

        See :meth:`OneCodexBase.where` for the full operator reference.
        """
        return super().where(
            *filters,
            sort=sort,
            limit=limit,
            public=public,
            filter=filter,
            updated_at=updated_at,
            date_collected=date_collected,
            date_sequenced=date_sequenced,
            description=description,
            external_sample_id=external_sample_id,
            library_type=library_type,
            location_lat=location_lat,
            location_lon=location_lon,
            location_string=location_string,
            name=name,
            platform=platform,
            sample_type=sample_type,
            starred=starred,
            sample=sample,
        )


class _SampleSchema(SampleSchema):
    _resource_path = "/api/v1/samples"
    _use_cursor_pagination = True

    # Use ApiRef for all reference fields
    owner: Union[Users, ApiRef]
    metadata: Union[Metadata, ApiRef]
    primary_classification: Optional[Union[Classifications, ApiRef]] = None
    project: Optional[Union[Projects, ApiRef]] = None
    tags: List[Union[Tags, ApiRef]] = []

    created_at: datetime
    updated_at: Optional[datetime] = None
    visibility: str

    @field_validator("created_at", "updated_at", mode="before")
    @classmethod
    def validate_datetime(cls, v):
        return datetime.fromisoformat(v)


class Samples(OneCodexBase, _SampleSchema, ResourceDownloadMixin):
    _allowed_methods = {
        "delete": None,
        "update": SampleUpdateSchema,
        "instances_public": None,
    }

    def __repr__(self):
        return '<{} {}: "{}">'.format(
            self.__class__.__name__, self.id, truncate_string(self.filename or "(N/A)", 24)
        )

    @classmethod
    def where(  # type: ignore[override]
        cls,
        *filters: str | dict,
        sort: str | list[str] | None = None,
        limit: int | None = None,
        public: bool = False,
        organization: bool = False,
        filter: Any = None,
        tags: list | None = None,
        created_at: datetime | DatetimeFilter = UNSET,
        updated_at: datetime | DatetimeFilter = UNSET,
        filename: str | StrFilter | None = UNSET,
        error_msg: str | StrFilter | None = UNSET,
        size: int | NumFilter | None = UNSET,
        status: str | StrFilter = UNSET,
        visibility: str | StrFilter = UNSET,
        metadata: Metadata | str | RefFilter = UNSET,
        owner: Users | str | RefFilter = UNSET,
        project: Projects | str | RefFilter | None = UNSET,
        **keyword_filters: Any,  # Metadata fields are forwarded; see body.
    ) -> "SampleCollection":
        """Query samples and return a :class:`SampleCollection`.

        Samples are the uploaded sequencing files (FASTA/FASTQ). Filter by
        any sample field, by metadata field (transparently joined), or by
        tags.

        Examples
        --------
        Find your recent FASTQ uploads::

            from datetime import datetime, timedelta, timezone
            since = (datetime.now(timezone.utc) - timedelta(days=7)).isoformat()
            ocx.Samples.where(
                created_at={"$gte": since},
                filename={"$iendswith": ".fastq.gz"},
            )

        Filter by tag (resolved by name, id, or :class:`Tags` instance)::

            ocx.Samples.where(tags=["trimmed", "human-depleted"])

        Filter by a metadata field — transparently joined::

            ocx.Samples.where(platform="Illumina NovaSeq 6000")

        Search public or organization-shared samples::

            ocx.Samples.where(filename={"$icontains": "hmp"}, public=True)
            ocx.Samples.where(organization=True)

        Apply a client-side predicate after fetching::

            ocx.Samples.where(filter=lambda s: not s.tags)

        Parameters
        ----------
        organization
            Search samples shared across your organization. Mutually
            exclusive with ``public``.
        public
            Search public samples (capped at 1000 results). Mutually
            exclusive with ``organization``.
        tags
            Tags to filter by. Accepts :class:`Tags` instances, tag ids, or
            tag names — all resolved to refs and combined with
            ``$containsall``.

        Returns
        -------
        :class:`SampleCollection`

        Notes
        -----
        Metadata ``custom`` fields aren't server-filterable — filter
        client-side via :meth:`SampleCollection.filter`.

        See :meth:`OneCodexBase.where` for the full operator reference.
        """
        from onecodex.models.collection import SampleCollection

        instances_route = "instances"
        if organization is True:
            instances_route = "instances_organization"
        if public is True:
            instances_route = "instances_public"

        effective_limit = limit if limit is not None else (1000 if public else None)

        # handle conversion of tag UUIDs or names to Tags objects
        if tags is None:
            tags = []
        if not isinstance(tags, list):
            tags = [tags]

        new_tags = []
        for t in tags:
            if isinstance(t, Tags):
                new_tags.append(t)
                continue

            if len(t) == 16:
                new_tag = Tags.get(t)
            else:
                where_tags = Tags.where(name=t)
                if len(where_tags) == 1:
                    new_tag = where_tags[0]
                elif len(where_tags) > 1:
                    raise OneCodexException("Multiple tags matched query: {}".format(t))
                else:
                    raise OneCodexException("Unknown tag specified: {}".format(t))
            new_tags.append(new_tag)

        # Merge inline field kwargs with leftover **keyword_filters.
        merged_filters = _drop_unset(
            created_at=created_at,
            updated_at=updated_at,
            filename=filename,
            error_msg=error_msg,
            size=size,
            status=status,
            visibility=visibility,
            metadata=metadata,
            owner=owner,
            project=project,
        )
        merged_filters.update(keyword_filters)
        if new_tags:
            merged_filters["tags"] = {"$containsall": new_tags}

        # we can only search metadata on our own samples currently
        # FIXME: we need to add `instances_public` and `instances_project` metadata routes to
        # mirror the ones on the samples
        md_search_keywords = {
            k: v
            for k, v in merged_filters.items()
            if k in Metadata.model_fields and k not in Samples.model_fields
        }
        if not public and not organization:
            if md_search_keywords:
                metadata_samples = [md.sample for md in Metadata.where(**md_search_keywords)]

        if md_search_keywords:
            # we tried searching by metadata fields
            if not metadata_samples:
                # there were no results, so don't bother with a slower query on Samples
                samples = []
            elif not (filters or merged_filters):
                # there were results, and there are no other filters to apply, so return them
                samples = metadata_samples
            else:
                # there were results, and we want to return the intersection of those with a query
                # on Samples using any non-metadata filters
                metadata_sample_ids = {s.id for s in metadata_samples}
                samples = super(Samples, cls).where(
                    *filters,
                    sort=sort,
                    limit=effective_limit,
                    filter=filter,
                    _instances=instances_route,
                    **merged_filters,
                )
                samples = [s for s in samples if s.id in metadata_sample_ids]
        else:
            # we did not try searching by metadata fields, so return whatever this gives us. in the
            # case that no filters/keyword_filters are specified, this is identical to Samples.all()
            samples = super(Samples, cls).where(
                *filters,
                sort=sort,
                limit=effective_limit,
                filter=filter,
                _instances=instances_route,
                **merged_filters,
            )

        return SampleCollection([s for s in samples[:effective_limit]], Samples)

    def save(self):
        """Send changes on this Samples object to the One Codex server.

        Changes to the metadata object and tags list are passed as well.
        """
        # TODO: Add a test for this
        if self.tags:
            for tag in self.tags:
                if tag.id is None:
                    tag.save()

        if self.project and not isinstance(self.project, Projects):
            try:
                self.project = get_project(self.project)
            except OneCodexException as e:
                raise OneCodexException("Error saving sample: {}".format(e))

        # FIXME: PATCH HERE
        super(Samples, self).save()

        if self.metadata is not None:
            self.metadata.save()

    @classmethod
    def preupload(cls, metadata=None, tags=None, project=None):
        """Create a sample in a waiting state where the files will be sent later on.

        Parameters
        ----------
        metadata : `dict`, optional

        tags : `list`, optional
            A list of optional tags to create. Tags must be passed as dictionaries with a single key
            `name` and the tag name, e.g., {"name": "my tag"}. New tags will be created on-the-fly
            as needed.

        project : `string`, optional
            UUID of project to associate this sample with.
        """
        project = get_project(project)

        resp = cls._client.post(
            f"{cls._api._base_url}{cls._resource_path}/preupload",
            json={
                "metadata": metadata,
                "tags": tags,
                "project": {"$ref": project.field_uri} if project else None,
            },
        )
        resp.raise_for_status()
        return resp.json()["sample_id"]

    @classmethod
    def upload(
        cls,
        files,
        metadata=None,
        tags=None,
        project=None,
        coerce_ascii=False,
        progressbar=None,
        sample_id=None,
        external_sample_id=None,
    ):
        """Upload a series of files to the One Codex server.

        Parameters
        ----------
        files : `string` or `tuple`
            A single path to a file on the system, or a tuple containing a pairs of paths. Tuple
            values  will be interleaved as paired-end reads and both files should contain the same
            number of records. Paths to single files will be uploaded as-is.

        metadata : `dict`, optional

        tags : `list`, optional
            A list of optional tags to create. Tags must be passed as dictionaries with a single key
            `name` and the tag name, e.g., {"name": "my tag"}. New tags will be created on-the-fly
            as needed.

        project : `string`, optional
            UUID of project to associate this sample with.

        coerce_ascii : `bool`, optional
            If true, rename unicode filenames to ASCII and issue warning.

        progressbar : `click.progressbar`, optional
            If passed, display a progress bar using Click.

        sample_id : `string`, optional
            If passed, will upload the file(s) to the sample with that id. Only works if the sample was pre-uploaded

        external_sample_id : `string`, optional
            If passed, will upload the file(s) to the sample with that metadata external id. Only works if the sample was pre-uploaded

        Returns
        -------
        A `Samples` object upon successful upload. None if the upload failed.
        """
        if not isinstance(files, str) and not isinstance(files, tuple):
            raise OneCodexException(
                "Please pass a string or tuple or forward and reverse filepaths."
            )

        if sample_id and external_sample_id:
            raise OneCodexException("Only pass sample_id OR external_sample_id, not both.")

        project = get_project(project)
        sample_id = upload_sequence(
            files,
            cls,
            metadata=metadata,
            tags=tags,
            project=project,
            coerce_ascii=coerce_ascii,
            progressbar=progressbar,
            sample_id=sample_id,
            external_sample_id=external_sample_id,
        )

        return cls.get(sample_id)
