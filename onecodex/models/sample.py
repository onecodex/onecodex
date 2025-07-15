from typing import Optional, List, Union
from datetime import datetime

from requests.exceptions import HTTPError

from onecodex.exceptions import OneCodexException
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


class MetadataSchema(MetadataSchema):
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


class Metadata(OneCodexBase, MetadataSchema):
    _resource_path = "/api/v1/metadata"
    _allowed_methods = {
        "update": MetadataPatchSchema,
    }


class _SampleSchema(SampleSchema):
    _resource_path = "/api/v1/samples"

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
    def where(cls, *filters, **keyword_filters):
        """Query and retrieve a set of samples.

        Parameters
        ----------
        limit : `int`, optional
            If set, retrieve a maximum of `limit` samples. Note a limit of 1000 samples is automatically
            enforced for queries of all public samples (`public=True`)
        organization : `bool`, optional
            If True, search all samples within your organization (including your samples). May not be combined with `public=True`.
        public : `bool`, optional
            If True, search all public samples (limited to 1000 results). May not be combined with `organization=True`.
        tags : `list`, optional
            A list of optional Tags to filter by. Tags should be `Tag` objects retrieved
            with `ocx.Tags.get()` or `ocx.Tags.where()`
        project : `Project`, optional
            Filter by a Project
        **keyword_filters : dict, optional
            Pass any additional sample or metadata attribute to filter by that attribute. Metadata filtering
            is *not* currently supported for

        Returns
        -------
        A `SampleCollection` object with samples matching the query
        """
        from onecodex.models.collection import SampleCollection

        public = keyword_filters.pop("public", False)
        organization = keyword_filters.pop("organization", False)
        instances_route = "instances"
        if organization is True:
            instances_route = "instances_organization"
        if public is True:
            instances_route = "instances_public"

        # Set default filters
        keyword_filters["_instances"] = instances_route
        keyword_filters.setdefault("limit", 1000 if public else None)

        # handle conversion of tag UUIDs or names to Tags objects
        tags = keyword_filters.pop("tags", [])

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

        if new_tags:
            keyword_filters["tags"] = {"$containsall": new_tags}

        # we can only search metadata on our own samples currently
        # FIXME: we need to add `instances_public` and `instances_project` metadata routes to
        # mirror the ones on the samples
        md_search_keywords = {
            k: v
            for k, v in keyword_filters.items()
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
            elif not (filters or keyword_filters):
                # there were results, and there are no other filters to apply, so return them
                samples = metadata_samples
            else:
                # there were results, and we want to return the intersection of those with a query
                # on Samples using any non-metadata filters
                metadata_sample_ids = {s.id for s in metadata_samples}
                samples = super(Samples, cls).where(*filters, **keyword_filters)
                samples = [s for s in samples if s.id in metadata_sample_ids]
        else:
            # we did not try searching by metadata fields, so return whatever this gives us. in the
            # case that no filters/keyword_filters are specified, this is identical to Samples.all()
            samples = super(Samples, cls).where(*filters, **keyword_filters)

        return SampleCollection([s for s in samples[: keyword_filters["limit"]]], Samples)

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
