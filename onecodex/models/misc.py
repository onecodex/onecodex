import warnings

from onecodex.exceptions import OneCodexException
from onecodex.lib.upload import upload_document
from onecodex.models import OneCodexBase
from onecodex.models.helpers import truncate_string, ResourceDownloadMixin


class Tags(OneCodexBase):
    _resource_path = "/api/v1/tags"

    def __init__(self, _resource=None, **kwargs):
        if kwargs:
            # try to lookup Tags with a where call using kwargs
            results = self.where(**kwargs)

            if len(results) == 0:
                super(Tags, self).__init__(_resource=_resource, **kwargs)
            elif len(results) == 1:
                self._resource = results[0]._resource
            elif len(results) > 1:
                raise OneCodexException("Multiple matches found for given criteria")
        else:
            super(Tags, self).__init__(_resource=_resource, **kwargs)

    def __repr__(self):
        return '<{} {}: "{}">'.format(
            self.__class__.__name__, self.id, truncate_string(self.name, 24)
        )

    def __hash__(self):
        return hash(self.name)


class Users(OneCodexBase):
    _resource_path = "/api/v1/users"


class Projects(OneCodexBase):
    _resource_path = "/api/v1/projects"

    @classmethod
    def search_public(cls, *filters, **keyword_filters):
        warnings.warn("Now supported via `.where(..., public=True)`", DeprecationWarning)
        keyword_filters["public"] = True
        keyword_filters["limit"] = 1000
        return cls.where(*filters, **keyword_filters)


class Jobs(OneCodexBase):
    _resource_path = "/api/v1/jobs"


class Documents(OneCodexBase, ResourceDownloadMixin):
    _resource_path = "/api/v1/documents"

    @classmethod
    def upload(cls, file_path, progressbar=None):
        """Uploads a series of files to the One Codex server.

        Parameters
        ----------
        file_path : `string`
            A path to a file on the system.
        progressbar : `click.progressbar`, optional
            If passed, display a progress bar using Click.

        Returns
        -------
        A `Samples` object upon successful upload. None if the upload failed.
        """
        res = cls._resource
        doc_id = upload_document(file_path, res._client.session, res, progressbar=progressbar)

        return cls.get(doc_id)
