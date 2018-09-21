import os
import sys
import warnings

import requests
from requests.exceptions import HTTPError
from six import string_types

from onecodex.exceptions import OneCodexException
from onecodex.models import OneCodexBase
from onecodex.models.helpers import truncate_string
from onecodex.lib.upload import upload  # upload_file


class Samples(OneCodexBase):
    _resource_path = '/api/v1/samples'

    def __repr__(self):
        return '<{} {}: "{}">'.format(self.__class__.__name__, self.id,
                                      truncate_string(self.filename, 24))

    @classmethod
    def where(cls, *filters, **keyword_filters):
        public = keyword_filters.get('public', False)
        instances_route = 'instances' if not public else 'instances_public'
        limit = keyword_filters.get('limit', None if not public else 1000)

        # we can only search metadata on our own samples currently
        # FIXME: we need to add `instances_public` and `instances_project` metadata routes to
        # mirror the ones on the samples
        metadata_samples = []
        if not public:
            md_schema = next(l for l in Metadata._resource._schema['links']
                             if l['rel'] == instances_route)

            md_where_schema = md_schema['schema']['properties']['where']['properties']
            md_search_keywords = {}
            for keyword in list(keyword_filters):
                # skip out on $uri to prevent duplicate field searches and the others to
                # simplify the checking below
                if keyword in ['$uri', 'sort', '_instances']:
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

        return samples[:limit]

    @classmethod
    def search_public(cls, *filters, **keyword_filters):
        warnings.warn('Now supported via `where(..., public=True)`', DeprecationWarning)
        keyword_filters['public'] = True
        keyword_filters['limit'] = 1000
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
    def upload(cls, filename, threads=None, validate=True, metadata=None, tags=None, project=None):
        """
        Uploads a series of files to the One Codex server. These files are automatically
        validated during upload.

        Parameters
        ----------
        path: list of strings or tuples
            List of full paths to the files. If one (or more) of the list items are a tuple, this
            is parsed as a set of files that are paired and the files are automatically
            iterleaved during upload.
        """
        # TODO: either raise/wrap UploadException or just us the new one in lib.samples
        # upload_file(filename, cls._resource._client.session, None, 100)
        res = cls._resource
        if isinstance(filename, string_types) or isinstance(filename, tuple):
            filename = [filename]
        samples = upload(filename, res._client.session, res, res._client._root_url + '/', threads=threads,
                         validate=validate, log_to=sys.stderr, metadata=metadata, tags=tags, project=project)
        return samples
        # FIXME: pass the auth into this so we can authenticate the callback?

    def download(self, path=None):
        """
        Downloads the original reads file (FASTA/FASTQ) from One Codex.

        Note that this may only work from within a notebook session and the file
        is not guaranteed to exist for all One Codex plan types.

        Parameters
        ----------
        path : string, optional
            Full path to save the file to. If omitted, defaults to the original filename
            in the current working directory.
        """
        if path is None:
            path = os.path.join(os.getcwd(), self.filename)
        try:
            url_data = self._resource.download_uri()
            resp = requests.get(url_data['download_uri'], stream=True)
            # TODO: use tqdm or ProgressBar here to display progress?
            with open(path, 'wb') as f_out:
                for data in resp.iter_content(chunk_size=1024):
                    f_out.write(data)
        except HTTPError as exc:
            if exc.response.status_code == 402:
                raise OneCodexException('You must either have a premium platform account or be in '
                                        'a notebook environment to download samples.')
            else:
                raise OneCodexException('Download failed with an HTTP status code {}.'.format(
                                        exc.response.status_code))


class Metadata(OneCodexBase):
    _resource_path = '/api/v1/metadata'

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
            if 'sample' in ref_props or '$uri' in ref_props:  # May not be there if not resolved!
                ref_props.pop('$uri', None)
                ref_props.pop('sample', None)
                self._resource._update(ref_props)
            else:
                super(Metadata, self).save()
