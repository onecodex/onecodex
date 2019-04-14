from requests.exceptions import HTTPError
from six import string_types
import warnings

from onecodex.exceptions import OneCodexException
from onecodex.lib.upload import upload_sequence
from onecodex.models import OneCodexBase, Projects, Tags
from onecodex.models.helpers import truncate_string, ResourceDownloadMixin


class Samples(OneCodexBase, ResourceDownloadMixin):
    _resource_path = "/api/v1/samples"

    def __repr__(self):
        return '<{} {}: "{}">'.format(
            self.__class__.__name__, self.id, truncate_string(self.filename, 24)
        )

    @classmethod
    def where(cls, *filters, **keyword_filters):
        from onecodex.models.collection import SampleCollection

        public = keyword_filters.get("public", False)
        instances_route = "instances" if not public else "instances_public"
        limit = keyword_filters.get("limit", None if not public else 1000)

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
        metadata_samples = []
        if not public:
            md_schema = next(
                l for l in Metadata._resource._schema["links"] if l["rel"] == instances_route
            )

            md_where_schema = md_schema["schema"]["properties"]["where"]["properties"]
            md_search_keywords = {}
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

            if len(md_search_keywords) > 0:
                metadata_samples = [md.sample for md in Metadata.where(**md_search_keywords)]

        samples = []
        if len(metadata_samples) == 0:
            samples = super(Samples, cls).where(*filters, **keyword_filters)

        if len(samples) > 0 and len(metadata_samples) > 0:
            # we need to filter samples to just include stuff from metadata_samples
            metadata_sample_ids = {s.id for s in metadata_samples}
            samples = [s for s in samples if s.id in metadata_sample_ids]
        elif len(metadata_samples) > 0:
            # we have to sort the metadata samples manually using the
            # sort parameters for the samples (and then the metadata parameters?)
            # TODO: implement this (see above block)
            samples = metadata_samples

        return SampleCollection([s._resource for s in samples[:limit]], Samples)

    @classmethod
    def search_public(cls, *filters, **keyword_filters):
        warnings.warn("Now supported via `where(..., public=True)`", DeprecationWarning)
        keyword_filters["public"] = True
        keyword_filters["limit"] = 1000
        return cls.where(*filters, **keyword_filters)

    def save(self):
        """
        Persist changes on this Samples object back to the One Codex server along with any changes
        on its metadata (if it has any).
        """
        super(Samples, self).save()
        if self.metadata is not None:
            self.metadata.save()

    @classmethod
    def upload(
        cls, files, metadata=None, tags=None, project=None, coerce_ascii=False, progressbar=None
    ):
        """Uploads a series of files to the One Codex server.

        Parameters
        ----------
        files : `string` or `tuple`
            A single path to a file on the system, or a tuple containing a pairs of paths. Tuple
            values  will be interleaved as paired-end reads and both files should contain the same
            number of records. Paths to single files will be uploaded as-is.
        metadata : `dict`, optional
        tags : `list`, optional
        project : `string`, optional
            UUID of project to associate this sample with.
        coerce_ascii : `bool`, optional
            If true, rename unicode filenames to ASCII and issue warning.
        progressbar : `click.progressbar`, optional
            If passed, display a progress bar using Click.

        Returns
        -------
        A `Samples` object upon successful upload. None if the upload failed.
        """
        res = cls._resource
        if not isinstance(files, string_types) and not isinstance(files, tuple):
            raise OneCodexException(
                "Please pass a string or tuple or forward and reverse filepaths."
            )

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
                raise OneCodexException("{} is not a valid project UUID".format(project))

            if isinstance(project_search, list):
                project = project_search[0]

        sample_id = upload_sequence(
            files,
            res._client.session,
            res,
            metadata=metadata,
            tags=tags,
            project=project,
            coerce_ascii=coerce_ascii,
            progressbar=progressbar,
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
