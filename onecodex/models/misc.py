import click
import os
import requests
import six
import warnings

from onecodex.exceptions import OneCodexException
from onecodex.lib.upload import upload_document
from onecodex.models import OneCodexBase
from onecodex.models.helpers import truncate_string


class Tags(OneCodexBase):
    _resource_path = '/api/v1/tags'

    def __init__(self, _resource=None, **kwargs):
        if kwargs:
            # try to lookup Tags with a where call using kwargs
            results = self.where(**kwargs)

            if len(results) == 0:
                super(Tags, self).__init__(_resource=_resource, **kwargs)
            elif len(results) == 1:
                self._resource = results[0]._resource
            elif len(results) > 1:
                raise OneCodexException('Multiple matches found for given criteria')
        else:
            super(Tags, self).__init__(_resource=_resource, **kwargs)

    def __repr__(self):
        return '<{} {}: "{}">'.format(self.__class__.__name__, self.id,
                                      truncate_string(self.name, 24))

    def __hash__(self):
        return hash(self.name)


class Users(OneCodexBase):
    _resource_path = '/api/v1/users'


class Projects(OneCodexBase):
    _resource_path = '/api/v1/projects'

    @classmethod
    def search_public(cls, *filters, **keyword_filters):
        warnings.warn('Now supported via `.where(..., public=True)`', DeprecationWarning)
        keyword_filters['public'] = True
        keyword_filters['limit'] = 1000
        return cls.where(*filters, **keyword_filters)


class Jobs(OneCodexBase):
    _resource_path = '/api/v1/jobs'


class Documents(OneCodexBase):
    _resource_path = '/api/v1/documents'

    def download(self, path=None, progressbar=False):
        """Downloads a document file from One Codex.

        Parameters
        ----------
        path : string, optional
            Full path to save the file to. If omitted, defaults to the original filename
            in the current working directory.
        """
        if path is None:
            path = os.path.join(os.getcwd(), self.filename)

        if os.path.exists(path):
            raise OneCodexException('{} already exists! Will not overwrite.'.format(path))

        try:
            url_data = self._resource.download_uri()
            resp = requests.get(url_data['download_uri'], stream=True)

            with open(path, 'wb') as f_out:
                if progressbar:
                    with click.progressbar(length=self.size, label=self.filename) as bar:
                        for data in resp.iter_content(chunk_size=1024):
                            bar.update(len(data))
                            f_out.write(data)
                else:
                    for data in resp.iter_content(chunk_size=1024):
                        f_out.write(data)
        except KeyboardInterrupt:
            os.remove(path)
        except requests.exceptions.HTTPError as exc:
            if exc.response.status_code == 401:
                raise OneCodexException('You must be logged in to download files.')
            elif exc.response.status_code == 402:
                raise OneCodexException('You must either have a premium platform account or be in '
                                        'a notebook environment to download files.')
            elif exc.response.status_code == 403:
                raise OneCodexException('You are not authorized to download this file.')
            else:
                raise OneCodexException('Download failed with an HTTP status code {}.'.format(
                                        exc.response.status_code))

        return path

    @classmethod
    def upload(cls, files, threads=None, log=None, progressbar=False):
        """Uploads a series of files to the One Codex server.

        Parameters
        ----------
        files : `string`, `tuple,` or `list`
            A list of paths to files on the system.
        threads : `integer`, optional
            Number of concurrent uploads. May provide a speedup.
        log : `logging.Logger`, optional
            Used to write status messages to a file or terminal.
        progressbar : `bool`, optional
            If true, display a progress bar using Click.
        """
        res = cls._resource
        if isinstance(files, six.string_types):
            files = [files]

        docs = upload_document(
            files, res._client.session, res, threads=threads, log=log, progressbar=progressbar
        )

        return docs
