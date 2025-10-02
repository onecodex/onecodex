import copy
import os
import os.path

import click
import requests

from onecodex.exceptions import OneCodexException


def coerce_search_value(search_value, field_name=None, base_class=None):
    from onecodex.models.base import OneCodexBase  # in here to prevent circular import

    if field_name == "$uri":
        return search_value.field_uri
    elif isinstance(search_value, OneCodexBase):
        return {"$ref": search_value.field_uri}
    elif isinstance(search_value, list):
        return [coerce_search_value(v) for v in search_value]
    elif isinstance(search_value, dict) and len(search_value) == 1:
        op = next(iter(search_value.keys()))
        value = search_value[op]
        return {op: coerce_search_value(value)}
    return search_value


def generate_potion_sort_clause(sort_items, sort_schema):
    if sort_items is None:
        return {}

    if not isinstance(sort_items, list):
        sort_items = [sort_items]

    sort_clause = {}
    for item in sort_items:
        if sort_schema and item.lstrip("^") not in sort_schema:
            raise AttributeError("Attribute {} cannot be sorted on".format(item.lstrip("^")))
        if item.startswith("^"):
            sort_clause[item[1:]] = False
        else:
            sort_clause[item] = True
    return sort_clause


def generate_potion_keyword_where(keyword_filters, where_schema, base_class):
    where = {}
    for keyword, search_value in keyword_filters.items():
        if keyword == "id":
            keyword = "$uri"

        if keyword not in base_class.model_fields:
            raise AttributeError("{} cannot be searched on {}".format(base_class.__name__, keyword))

        if isinstance(search_value, list):
            where[keyword] = {
                "$in": [coerce_search_value(v, keyword, base_class) for v in search_value]
            }
        else:
            where[keyword] = coerce_search_value(search_value, keyword, base_class)

    return where


def truncate_string(s, length=24):
    if len(s) < length - 3:
        return s
    else:
        s = s[0 : (length - 3)]
        if s[-1] == ".":
            s = s + ".."
        else:
            s = s + "..."
        return s


class ResourceDownloadMixin(object):
    def download(self, path=None, file_obj=None, progressbar=False):
        """Download files from One Codex.

        Parameters
        ----------
        path : `string`, optional
            Full path to save the file to. If omitted, defaults to the original filename
            in the current working directory.
        file_obj : file-like object, optional
            Rather than save the file to a path, write it to this file-like object.
        progressbar : `bool`
            Display a progress bar using Click for the download?

        Returns
        -------
        `string`
            The path the file was downloaded to, if applicable. Otherwise, None.

        Notes
        -----
        If no arguments specified, defaults to download the file as the original filename
        in the current working directory. If `file_obj` given, will write data into the
        passed file-like object. If `path` given, will download the file to the path provided,
        but will not overwrite any existing files.
        """
        return self._download(
            "download_uri",
            _filename=self.filename,
            path=path,
            file_obj=file_obj,
            progressbar=progressbar,
        )

    def _download(
        self,
        _resource_method,
        _filename=None,
        use_client_session=False,
        path=None,
        file_obj=None,
        progressbar=False,
    ):
        from requests.adapters import HTTPAdapter
        from requests.packages.urllib3.util.retry import Retry

        if hasattr(self, "visibility") and self.visibility == "awaiting data":
            raise OneCodexException("Sample has not finished processing. Please try again later.")

        if path and file_obj:
            raise OneCodexException("Please specify only one of: path, file_obj")

        try:
            if path is None and file_obj is None and _filename is None:
                raise OneCodexException("Please specify `path`, `file_obj`, or `_filename`.")

            if path is None and file_obj is None:
                path = os.path.join(os.getcwd(), _filename)

            if path and os.path.exists(path):
                raise OneCodexException("{} already exists. Will not overwrite.".format(path))

            download_link_info = self._client.post(
                f"{self._api._base_url}{self.field_uri}/{_resource_method}"
            )
            download_link_info.raise_for_status()
            link = download_link_info.json()["download_uri"]

            if use_client_session:
                session = copy.deepcopy(self._client.session)
            else:
                # Create a new session with custom retries for downloads. Do not use HTTPClient
                # because this request goes to S3 and uses a different auth mechanism
                session = requests.Session()

            # Retry up to 5 times with backoff timing of 2s, 4s, 8s, 16s, and 32s (applies to all
            # HTTP methods). 404 is included for cases where the file is being asynchronously
            # uploaded to S3 and is expected to be available soon.
            retry_strategy = Retry(
                total=5,
                backoff_factor=2,
                status_forcelist=[404, 429, 500, 502, 503, 504],
                allowed_methods=None,
            )
            adapter = HTTPAdapter(max_retries=retry_strategy)

            session.mount("http://", adapter)
            session.mount("https://", adapter)
            resp = session.get(link, stream=True)

            if path:
                f_out = open(path, "wb")
            else:
                f_out = file_obj

            if progressbar:
                progress_label = os.path.basename(path) if path else self.filename
                with click.progressbar(length=self.size, label=progress_label) as bar:
                    for data in resp.iter_content(chunk_size=1024):
                        bar.update(len(data))
                        f_out.write(data)
            else:
                for data in resp.iter_content(chunk_size=1024):
                    f_out.write(data)

            # do not close the handle if file_obj is used
            if not file_obj:
                f_out.close()

        except KeyboardInterrupt:
            if path and os.path.exists(path):
                os.remove(path)
            raise
        except requests.exceptions.HTTPError as exc:
            if exc.response.status_code == 401:
                raise OneCodexException("You must be logged in to download files.")
            elif exc.response.status_code == 402:
                raise OneCodexException(
                    "You must either have a premium platform account or be in "
                    "a notebook environment to download files. Please feel free to contact us "
                    "about your subscription at support@onecodex.com."
                )
            elif exc.response.status_code == 403:
                raise OneCodexException("You are not authorized to download this file.")
            else:
                raise OneCodexException(
                    "Download failed with an HTTP status code {}.".format(exc.response.status_code)
                )

        return path
