from requests.exceptions import HTTPError
from six import string_types
import warnings

from onecodex.exceptions import OneCodexException
from onecodex.lib.upload import upload_sequence, preupload_sample
from onecodex.models import OneCodexBase, Projects, Tags
from onecodex.models.helpers import truncate_string, ResourceDownloadMixin


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


class Samples(OneCodexBase, ResourceDownloadMixin):
    _resource_path = "/api/v1/samples"

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
            else:
                # is it a uuid?
                if len(t) == 16:
                    tag_obj = Tags.get(t)

                    if tag_obj is not None:
                        new_tags.append(tag_obj)
                        continue

                # is it a name?
                tag_obj = Tags.where(name=t)

                if len(tag_obj) == 1:
                    new_tags.append(tag_obj[0])
                    continue
                elif len(tag_obj) > 1:
                    raise OneCodexException("Multiple tags matched query: {}".format(t))

                raise OneCodexException("Unknown tag specified: {}".format(t))

        if new_tags:
            keyword_filters["tags"] = new_tags

        # we can only search metadata on our own samples currently
        # FIXME: we need to add `instances_public` and `instances_project` metadata routes to
        # mirror the ones on the samples
        md_search_keywords = {}
        if not public and not organization:
            md_schema = next(
                link
                for link in Metadata._resource._schema["links"]
                if link["rel"] == instances_route
            )

            md_where_schema = md_schema["schema"]["properties"]["where"]["properties"]

            for keyword in list(keyword_filters):
                # skip out on $uri to prevent duplicate field searches and the others to
                # simplify the checking below
                if keyword in ["$uri", "sort", "_instances"]:
                    continue
                elif keyword in md_where_schema:
                    md_search_keywords[keyword] = keyword_filters.pop(keyword)

            # TODO: should one be able to sort on metadata? here and on the merged list?
            # md_sort_schema = md_schema['schema']['properties']['sort']['properties']
            # # pull out any metadata sort parameters
            # sort = keyword_filters.get('sort', [])
            # if not isinstance(sort, list):
            #     sort = [sort]
            # passthrough_sort = []
            # for keyword in sort:
            #     if keyword in md_sort_schema:
            #         # TODO: set up sort for metadata
            #         pass
            #     else:
            #         passthrough_sort.append(keyword)
            # keyword_filters['sort'] = passthrough_sort

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

        return SampleCollection([s._resource for s in samples[: keyword_filters["limit"]]], Samples)

    @classmethod
    def search_public(cls, *filters, **keyword_filters):
        warnings.warn("Now supported via `where(..., public=True)`", DeprecationWarning)
        keyword_filters["public"] = True
        return cls.where(*filters, **keyword_filters)

    def save(self):
        """Send changes on this Samples object to the One Codex server.

        Changes to the metadata object and tags list are passed as well.
        """
        # any newly-created tags should be associated with this sample and saved
        if self.tags:
            for tag in self.tags:
                if tag.id is None:
                    tag._resource.__dict__["_Reference__properties"]["sample"] = self._resource
                    tag._resource.save()

        if self.project and not isinstance(self.project, Projects):
            try:
                self.project = get_project(self.project)
            except OneCodexException as e:
                raise OneCodexException("Error saving sample: {}".format(e))

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
        res = cls._resource
        project = get_project(project)

        sample_id = preupload_sample(res, metadata, tags, project)
        return cls.get(sample_id)

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
        if not isinstance(files, string_types) and not isinstance(files, tuple):
            raise OneCodexException(
                "Please pass a string or tuple or forward and reverse filepaths."
            )

        if sample_id and external_sample_id:
            raise OneCodexException("Only pass sample_id OR external_sample_id, not both.")

        project = get_project(project)

        sample_id = upload_sequence(
            files,
            cls._resource,
            metadata=metadata,
            tags=tags,
            project=project,
            coerce_ascii=coerce_ascii,
            progressbar=progressbar,
            sample_id=sample_id,
            external_sample_id=external_sample_id,
        )

        return cls.get(sample_id)

    def __hash__(self):
        return hash(self.id)


class Metadata(OneCodexBase):
    _resource_path = "/api/v1/metadata"

    def save(self):
        if self.id is None:
            super(Metadata, self).save()  # Create
        else:  # Update
            # Hack: Sample is read and create-only
            # but Potion will try to update since it's not marked
            # readOnly in the schema; we also make sure
            # the linked metadata object is resolved since
            # we auto-save it alongside the sample
            if self._resource._uri and self._resource._status is None:
                assert isinstance(self._resource._properties, dict)

            # Then eject samplea and uri as needed
            ref_props = self._resource._properties
            if "sample" in ref_props or "$uri" in ref_props:  # May not be there if not resolved!
                ref_props.pop("$uri", None)
                ref_props.pop("sample", None)
                self._resource._update(ref_props)
            else:
                super(Metadata, self).save()
