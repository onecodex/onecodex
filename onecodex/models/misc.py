import warnings
from typing import Optional, Union

from onecodex.exceptions import MethodNotSupported
from onecodex.lib.upload import upload_document

# from onecodex.lib.upload import upload_asset
from onecodex.models.helpers import truncate_string, ResourceDownloadMixin

from onecodex.models.base import OneCodexBase, ApiRef
from onecodex.models.generated import TagSchema as GeneratedTagSchema
from onecodex.models.generated import UserSchema as GeneratedUserSchema
from onecodex.models.generated import ProjectSchema as GeneratedProjectSchema
from onecodex.models.generated import JobSchema as GeneratedJobSchema
from onecodex.models.generated import DocumentSchema as GeneratedDocumentSchema
# from onecodex.models.generated import AssetSchema as GeneratedAssetSchema


class Tags(OneCodexBase, GeneratedTagSchema):
    _resource_path = "/api/v1/tags"

    # TODO: Implement this logic
    # def __init__(self, _resource=None, name=None, sample=None):
    #     if name:
    #         # try to lookup Tags with a where call using kwargs
    #         results = self.where(name=name)

    #         if len(results) == 0:
    #             super(Tags, self).__init__(name=name, sample=sample)
    #         elif len(results) == 1:
    #             # FIXME: Remove _resource reference
    #             self._resource = results[0]._resource
    #         elif len(results) > 1:
    #             raise OneCodexException("Multiple matches found for given criteria")
    #     else:
    #         # FIXME: Remove _resource reference
    #         super(Tags, self).__init__(_resource=_resource)

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


class Users(OneCodexBase, GeneratedUserSchema):
    _resource_path = "/api/v1/users"


class Projects(OneCodexBase, GeneratedProjectSchema):
    _resource_path = "/api/v1/projects"

    owner: Union[Users, ApiRef]

    def delete(self):
        self._client.delete(f"{self._api._base_url}{self.field_uri}")

    @classmethod
    def search_public(cls, *filters, **keyword_filters):
        warnings.warn("Now supported via `.where(..., public=True)`", DeprecationWarning)
        keyword_filters["public"] = True
        keyword_filters["limit"] = 1000
        return cls.where(*filters, **keyword_filters)


class Jobs(OneCodexBase, GeneratedJobSchema):
    _resource_path = "/api/v1/jobs"


class Documents(OneCodexBase, GeneratedDocumentSchema, ResourceDownloadMixin):
    _resource_path = "/api/v1/documents"
    size: Optional[int] = None

    def delete(self):
        self._client.delete(f"{self._api._base_url}{self.field_uri}")

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
        doc_id = upload_document(file_path, cls._resource, progressbar=progressbar)

        return cls.get(doc_id)


# Not supported in OpenAPI yet...
# class Assets(OneCodexBase, GeneratedAssetSchema, ResourceDownloadMixin):
#     _resource_path = "/api/v1_experimental/assets"

#     @classmethod
#     def upload(cls, file_path, progressbar=None, name=None):
#         """Upload a file to an asset.

#         Parameters
#         ----------
#         file_path : `string`
#             A path to a file on the system.
#         progressbar : `click.progressbar`, optional
#             If passed, display a progress bar using Click.
#         name : `string`, optional
#             If passed, name is sent with upload request and is associated with asset.

#         Returns
#         -------
#         An `Assets` object upon successful upload. None if the upload failed.
#         """
#         asset_id = upload_asset(file_path, cls._resource, progressbar=progressbar, name=name)

#         return cls.get(asset_id)
