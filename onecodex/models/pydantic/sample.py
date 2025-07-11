from typing import Optional, List, get_origin, get_args, Union
from datetime import datetime

from onecodex.models.pydantic.base import ApiBaseModel
from onecodex.models.pydantic.analysis import Classifications

from onecodex.models.helpers import truncate_string
from onecodex.models.generated import SampleSchema as GeneratedSampleSchema


from pydantic import field_validator

from onecodex.models.generated import MetadataSchema as GeneratedMetadataSchema
from onecodex.models.pydantic.misc import Users, Projects, Tags


class MetadataSchema(GeneratedMetadataSchema):
    @field_validator("platform", "library_type", "sample_type")
    @classmethod
    def convert_enum_to_str(cls, v):
        if v is None:
            return v
        return str(v.value) if hasattr(v, "value") else v


class Metadata(ApiBaseModel, MetadataSchema):
    _resource_path = "/api/v1/metadata"

    def save(self):
        # TODO: Implement any custom logic required for metadata
        raise NotImplementedError

    #     if self.id is None:
    #         super(Metadata, self).save()  # Create
    #     else:  # Update
    #         # Hack: Sample is read and create-only
    #         # but Potion will try to update since it's not marked
    #         # readOnly in the schema; we also make sure
    #         # the linked metadata object is resolved since
    #         # we auto-save it alongside the sample
    #         if self._resource._uri and self._resource._status is None:
    #             assert isinstance(self._resource._properties, dict)

    #         # Then eject samplea and uri as needed
    #         ref_props = self._resource._properties
    #         if "sample" in ref_props or "$uri" in ref_props:  # May not be there if not resolved!
    #             ref_props.pop("$uri", None)
    #             ref_props.pop("sample", None)
    #             self._resource._update(**ref_props)
    #         else:
    #             super(Metadata, self).save()


class _SampleSchema(GeneratedSampleSchema):
    _resource_path = "/api/v1/samples"

    # TODO: Remap to our classes here vs. Pydantic generated ones
    owner: Users
    metadata: Metadata
    primary_classification: Optional[Classifications] = None
    project: Optional[Projects] = None
    tags: List[Tags] = []

    created_at: datetime
    updated_at: Optional[datetime] = None

    @classmethod
    def _resolve_single_ref(cls, value, field_annotation):
        """Resolve a single $ref dict to an object via a GET request."""
        if isinstance(value, dict) and len(value) == 1 and "$ref" in value:
            ref_uri = value["$ref"]

            # Prevent circular references: don't resolve references back to samples
            if cls._resource_path in ref_uri:
                return value  # Keep as raw $ref dict

            # Handle Optional types (Union[SomeClass, None])
            if get_origin(field_annotation) is Union:
                # Get the non-None type from Union
                args = get_args(field_annotation)
                target_class = next(arg for arg in args if arg is not type(None))
            elif get_origin(field_annotation) is list:
                # Handle List[SomeClass] - get the inner type
                target_class = get_args(field_annotation)[0]
            else:
                target_class = field_annotation

            ref_id = ref_uri.split("/")[-1]
            return target_class.get(ref_id)
        return value

    @field_validator(
        "owner", "metadata", "primary_classification", "project", "tags", mode="before"
    )
    @classmethod
    def _resolve_ref(cls, v, info):
        """Resolve $ref fields automatically to objects via a GET request. Prevents circular references to samples."""
        # Get the field annotation
        field_name = info.field_name
        field_annotation = cls.model_fields[field_name].annotation

        if isinstance(v, list):
            # Handle arrays (like tags)
            return [cls._resolve_single_ref(item, field_annotation) for item in v]
        else:
            # Handle single values
            return cls._resolve_single_ref(v, field_annotation)

    @field_validator("created_at", "updated_at", mode="before")
    @classmethod
    def validate_datetime(cls, v):
        return datetime.fromisoformat(v)


class Samples(ApiBaseModel, _SampleSchema):
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
        raise NotImplementedError

    def save(self):
        """Send changes on this Samples object to the One Codex server.

        Changes to the metadata object and tags list are passed as well.
        """
        raise NotImplementedError

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
        raise NotImplementedError

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
        raise NotImplementedError
