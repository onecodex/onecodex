import click
import inspect
import os
import os.path
import requests

from onecodex.exceptions import OneCodexException, UnboundObject


def as_uri(uuid, base_class):
    return base_class._resource._schema._uri.rstrip("#") + "/" + uuid


def coerce_search_value(search_value, field_name, base_class):
    from onecodex.models import OneCodexBase  # in here to prevent circular import

    if field_name == "$uri":
        return as_uri(field_name, base_class)
    elif isinstance(search_value, OneCodexBase):
        return {"$ref": search_value._resource._uri}
    return search_value


def check_bind(self_or_cls):
    if not hasattr(self_or_cls, "_resource"):
        name = "class" if inspect.isclass(self_or_cls) else "instance"
        raise UnboundObject("This {} is not associated with an API binding.".format(name))


def generate_potion_sort_clause(sort_items, sort_schema):
    if sort_items is None:
        return {}

    if not isinstance(sort_items, list):
        sort_items = [sort_items]

    sort_clause = {}
    for item in sort_items:
        if item.lstrip("^") not in sort_schema:
            raise AttributeError("Attribute {} can not be sorted on".format(item.lstrip("^")))
        if item.startswith("^"):
            sort_clause[item[1:]] = False
        else:
            sort_clause[item] = True
    return sort_clause


def generate_potion_keyword_where(keyword_filters, where_schema, base_class):
    where = {}
    for keyword in keyword_filters:
        search_value = keyword_filters[keyword]

        if keyword == "id":
            keyword = "$uri"

        if keyword not in where_schema:
            raise AttributeError(
                "{} can not be searched on {}".format(base_class.__name__, keyword)
            )

        avail_searches = [v["required"] for v in where_schema[keyword]["anyOf"] if "required" in v]
        # flatten the list
        avail_searches = [item for sublist in avail_searches for item in sublist]

        # TODO: do schema type checking here too?
        if "$eq" not in avail_searches and "$containsall" in avail_searches:
            if not isinstance(search_value, list):
                search_value = [search_value]
            where[keyword] = {
                "$containsall": [coerce_search_value(v, keyword, base_class) for v in search_value]
            }
        elif isinstance(search_value, list):
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
            self.filename,
            use_potion_session=False,
            path=path,
            file_obj=file_obj,
            progressbar=progressbar,
        )

    def _download(
        self,
        _resource_method,
        _filename=None,
        use_potion_session=False,
        path=None,
        file_obj=None,
        progressbar=False,
    ):
        from requests.adapters import HTTPAdapter
        from requests.packages.urllib3.util.retry import Retry

        if hasattr(self._resource, "visibility") and self._resource.visibility == "awaiting data":
            raise OneCodexException("Sample has not finished processing. Please try again later.")

        if path and file_obj:
            raise OneCodexException("Please specify only one of: path, file_obj")

        try:
            method_to_call = getattr(self._resource, _resource_method)
            download_link_info = method_to_call()

            if path is None and file_obj is None:
                if _filename is None:
                    if "save_as_filename" not in download_link_info:
                        raise OneCodexException(
                            "Please specify `path`, `file_obj`, or `_filename`."
                        )
                    _filename = download_link_info["save_as_filename"]
                path = os.path.join(os.getcwd(), _filename)

            if path and os.path.exists(path):
                raise OneCodexException("{} already exists. Will not overwrite.".format(path))

            if use_potion_session:
                session = self._resource._client.session
            else:
                session = requests.Session()

            link = download_link_info["download_uri"]

            # Retry up to 5 times with backoff timing of 2s, 4s, 8s, 16s, and 32s (applies to all
            # HTTP methods). 404 is included for cases where the file is being asynchronously
            # uploaded to S3 and is expected to be available soon.
            retry_strategy = Retry(
                total=5,
                backoff_factor=2,
                status_forcelist=[404, 429, 500, 502, 503, 504],
                method_whitelist=False,
            )
            adapter = HTTPAdapter(max_retries=retry_strategy)
            session.mount("http://", adapter)
            session.mount("https://", adapter)

            resp = session.get(link, stream=True)

            with (open(path, "wb") if path else file_obj) as f_out:
                if progressbar:
                    progress_label = os.path.basename(path) if path else self.filename
                    with click.progressbar(length=self.size, label=progress_label) as bar:
                        for data in resp.iter_content(chunk_size=1024):
                            bar.update(len(data))
                            f_out.write(data)
                else:
                    for data in resp.iter_content(chunk_size=1024):
                        f_out.write(data)
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
