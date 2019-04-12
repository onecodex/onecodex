from __future__ import print_function, division

import bz2
from collections import deque, OrderedDict
import gzip
import logging
import os
import re
import requests
from requests_toolbelt import MultipartEncoder
import six
import time
from unidecode import unidecode
import warnings

from onecodex.exceptions import (
    OneCodexException,
    UploadException,
    RetryableUploadException,
    raise_connectivity_error,
    raise_api_error,
)
from onecodex.utils import atexit_register, atexit_unregister, snake_case


DEFAULT_THREADS = 4


# buffer code from
# http://stackoverflow.com/questions/2192529/python-creating-a-streaming-gzipd-file-like/2193508
class Buffer(object):
    def __init__(self):
        self._buf = deque()
        self._size = 0
        self.closed = False

    def __len__(self):
        return self._size

    def write(self, data):
        self._buf.append(data)
        self._size += len(data)

    def read(self, size=-1):
        if size < 0:
            size = self._size
        ret_list = []
        while size > 0 and len(self._buf):
            s = self._buf.popleft()
            size -= len(s)
            ret_list.append(s)
        if size < 0:
            ret_list[-1], remainder = ret_list[-1][:size], ret_list[-1][size:]
            self._buf.appendleft(remainder)
        ret = b"".join(ret_list)
        self._size -= len(ret)
        return ret

    def flush(self):
        pass

    def close(self):
        self.closed = True


class FASTXInterleave(object):
    """Wrapper around two `file` objects that decompresses gzip or bz2, where applicable, and
    interleaves the two files either two or four lines at a time Yields uncompressed data.

    Parameters
    ----------
    file_path : `string`
        Path to file.
    file_size : `integer`
        Accurate size of file on disk. Used by `requests_toolbelt.MultipartEncoder` when
        determining how much data is left to read.
    file_format : {'fasta', 'fastq'}, optional
        Determines how many lines to read from each file. FASTA reads two, FASTQ reads four.
    progressbar : `click.progressbar`, optional
        The progress bar to update.
    """

    def __init__(self, file_path, file_size, file_format="fastq", progressbar=None):
        if file_path[0].endswith(".gz") or file_path[1].endswith(".gz"):
            self._fp_left = gzip.GzipFile(file_path[0], mode="rb")
            self._fp_right = gzip.GzipFile(file_path[1], mode="rb")
        elif file_path[0].endswith(".bz2") or file_path[1].endswith(".bz2"):
            self._fp_left = bz2.BZ2File(file_path[0], mode="rb")
            self._fp_right = bz2.BZ2File(file_path[1], mode="rb")
        else:
            self._fp_left = open(file_path[0], mode="rb")
            self._fp_right = open(file_path[1], mode="rb")

        if file_format == "fasta":
            raise OneCodexException("Interleaving FASTA files is currently unsupported")
        elif file_format == "fastq":
            self._lines_per_record = 4
        else:
            raise OneCodexException("file_format must be one of: fastq, fasta")

        self._tell = 0
        self._fsize = file_size
        self._buf = Buffer()

        self.progressbar = progressbar
        self.mime_type = "text/plain"

    @property
    def len(self):
        """Size of data left to be read.

        Notes
        -----
        `_fsize` is the uncompressed size of both files in the pair, summed. `_tell` is the amount
        of uncompressed data that has been read from this wrapper thus far.
        """
        return self._fsize - self._tell

    def read(self, size=-1):
        while len(self._buf) < size or size < 0:
            for fp in [self._fp_left, self._fp_right]:
                count = 0
                for line in fp:
                    self._buf.write(line)
                    count += 1
                    if count == self._lines_per_record:
                        break

            if count == 0:
                break

        bytes_read = self._buf.read(size)
        self._tell += len(bytes_read)

        if self.progressbar:
            self.progressbar.update(len(bytes_read))

        return bytes_read

    def seek(self, loc):
        """Called if upload fails and must be retried."""
        assert loc == 0

        # rewind progress bar
        if self.progressbar:
            self.progressbar.update(-self._tell)

        self._fp_left.seek(loc)
        self._fp_right.seek(loc)
        self._tell = loc
        self._buf = Buffer()

    def close(self):
        self._fp_left.close()
        self._fp_right.close()


class FilePassthru(object):
    """Wrapper around `file` object that updates a progress bar and guesses mime-type.

    Parameters
    ----------
    file_path : `string`
        Path to file.
    file_size : `integer`
        Accurate size of file on disk. Used by `requests_toolbelt.MultipartEncoder` when
        determining how much data is left to read.
    progressbar : `click.progressbar`, optional
        The progress bar to update.
    """

    def __init__(self, file_path, file_size, progressbar=None):
        self._fp = open(file_path, mode="rb")
        self._fsize = file_size
        self.progressbar = progressbar

        _, ext = os.path.splitext(file_path)

        if ext in {".gz", ".gzip"}:
            self.mime_type = "application/x-gzip"
        elif ext in {".bz", ".bz2", ".bzip", ".bzip2"}:
            self.mime_type = "application/x-bzip2"
        else:
            self.mime_type = "text/plain"

    def read(self, size=-1):
        bytes_read = self._fp.read(size)

        if self.progressbar:
            self.progressbar.update(len(bytes_read))

        return bytes_read

    @property
    def len(self):
        """Size of data left to be read."""
        return self._fsize - self._fp.tell()

    def seek(self, loc):
        """Called if upload fails and must be retried."""
        assert loc == 0

        # rewind progress bar
        if self.progressbar:
            self.progressbar.update(-self._fp.tell())

        self._fp.seek(loc)

    def close(self):
        self._fp.close()


def interleaved_filename(file_path):
    """Return filename used to represent a set of paired-end files. Assumes Illumina-style naming
    conventions where each file has _R1_ or _R2_ in its name."""
    if not isinstance(file_path, tuple):
        raise OneCodexException("Cannot get the interleaved filename without a tuple.")
    if re.match(".*[._][Rr][12][_.].*", file_path[0]):
        return re.sub("[._][Rr][12]", "", file_path[0])
    else:
        warnings.warn("Paired-end filenames do not match--are you sure they are correct?")
        return file_path[0]


def _file_size(file_path, uncompressed=False):
    """Return size of a single file, compressed or uncompressed"""
    _, ext = os.path.splitext(file_path)

    if uncompressed:
        if ext in {".gz", ".gzip"}:
            with gzip.GzipFile(file_path, mode="rb") as fp:
                try:
                    fp.seek(0, os.SEEK_END)
                    return fp.tell()
                except ValueError:
                    # on python2, cannot seek from end and must instead read to end
                    fp.seek(0)
                    while len(fp.read(8192)) != 0:
                        pass
                    return fp.tell()
        elif ext in {".bz", ".bz2", ".bzip", ".bzip2"}:
            with bz2.BZ2File(file_path, mode="rb") as fp:
                fp.seek(0, os.SEEK_END)
                return fp.tell()

    return os.path.getsize(file_path)


def _file_stats(file_path, enforce_fastx=True):
    """Return information about the file path (or paths, if paired), prior to upload.

    Parameters
    ----------
    file_path : `string` or `tuple`
        System path to the file(s) to be uploaded

    Returns
    -------
    `string`
        Filename, minus compressed extension (.gz or .bz2). If paired, use first path to generate
        the filename that will be used to represent both paths in the pair.
    `integer`
        If paired, the uncompressed file size of both files in the path. If single, the raw file
        size whether compressed or not. Pairs are always uploaded uncompressed, whereas single files
        are uploaded in whatever format they're in. One Codex will uncompress and re-compress as
        appropriate.
    {'fasta', 'fastq'}
        The format of the file being uploaded, guessed only by its extension. If paired, this
        determines how many lines to pull from each file during interleaving.
    """
    if isinstance(file_path, tuple):
        assert len(file_path) == 2
        file_size = sum(_file_size(f, uncompressed=True) for f in file_path)
        file_path = interleaved_filename(file_path)
        paired = True
    else:
        file_size = _file_size(file_path, uncompressed=False)
        paired = False

    new_filename, ext = os.path.splitext(os.path.basename(file_path))

    if ext in {".gz", ".gzip", ".bz", ".bz2", ".bzip"}:
        compressed = ext
        new_filename, ext = os.path.splitext(new_filename)
    else:
        compressed = ""

    # strip compressed extension if paired-end, since we're going to upload uncompressed
    if paired and compressed:
        final_filename = new_filename + ext
    else:
        final_filename = new_filename + ext + compressed

    if enforce_fastx:
        if ext in {".fa", ".fna", ".fasta"}:
            file_format = "fasta"
        elif ext in {".fq", ".fastq"}:
            file_format = "fastq"
        else:
            raise UploadException(
                "{}: extension must be one of .fa, .fna, .fasta, .fq, .fastq".format(final_filename)
            )
    else:
        file_format = None

    if file_size == 0:
        raise UploadException("{}: empty files can not be uploaded".format(final_filename))

    return final_filename, file_size, file_format


# this lets us turn off the click progressbar context manager and is python2 compatible
# https://stackoverflow.com/questions/45187286/how-do-i-write-a-null-no-op-contextmanager-in-python
class FakeProgressBar(object):
    pct = 0
    label = ""

    def __init__(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def finish(self):
        pass

    def update(self, size):
        pass


def _call_init_upload(file_name, file_size, metadata, tags, project, samples_resource):
    """Call init_upload at the One Codex API and return data used to upload the file.

    Parameters
    ----------
    file_name : `string`
        The file_name you wish to associate this fastx file with at One Codex.
    file_size : `integer`
        Accurate size of file to be uploaded, in bytes.
    metadata : `dict`, optional
    tags : `list`, optional
    project : `string`, optional
        UUID of project to associate this sample with.
    samples_resource : `onecodex.models.Samples`
        Wrapped potion-client object exposing `init_upload` and `confirm_upload` routes to mainline.

    Returns
    -------
    `dict`
        Contains, at a minimum, 'upload_url' and 'sample_id'. Should also contain various additional
        data used to upload the file to fastx-proxy, a user's S3 bucket, or an intermediate bucket.
    """
    upload_args = {
        "filename": file_name,
        "size": file_size,
        "upload_type": "standard",  # this is multipart form data
    }

    if metadata:
        # format metadata keys as snake case
        new_metadata = {}

        for md_key, md_val in metadata.items():
            new_metadata[snake_case(md_key)] = md_val

        upload_args["metadata"] = new_metadata

    if tags:
        upload_args["tags"] = tags

    if project:
        upload_args["project"] = getattr(project, "id", project)

    try:
        upload_info = samples_resource.init_upload(upload_args)
    except requests.exceptions.HTTPError as e:
        raise_api_error(e.response, state="init")
    except requests.exceptions.ConnectionError:
        raise_connectivity_error(file_name)

    return upload_info


def _make_retry_fields(file_name, metadata, tags, project):
    """Generate fields to send to init_multipart_upload in the case that a Sample upload via
    fastx-proxy fails.

    Parameters
    ----------
    file_name : `string`
        The file_name you wish to associate this fastx file with at One Codex.
    metadata : `dict`, optional
    tags : `list`, optional
    project : `string`, optional
        UUID of project to associate this sample with.

    Returns
    -------
    `dict`
        Contains metadata fields that will be integrated into the Sample model created when
        init_multipart_upload is called.
    """
    upload_args = {"filename": file_name}

    if metadata:
        # format metadata keys as snake case
        new_metadata = {}

        for md_key, md_val in metadata.items():
            new_metadata[snake_case(md_key)] = md_val

        upload_args["metadata"] = new_metadata

    if tags:
        upload_args["tags"] = tags

    if project:
        upload_args["project"] = getattr(project, "id", project)

    return upload_args


def upload_sequence(
    files,
    session,
    samples_resource,
    metadata=None,
    tags=None,
    project=None,
    coerce_ascii=False,
    progressbar=None,
):
    """Uploads a sequence file (or pair of files) to the One Codex server via either our proxy or directly to S3.

    Parameters
    ----------
    files : `list`
        A list of paths to files on the system, or tuples containing pairs of paths. Tuples will be
        interleaved as paired-end reads and both files should contain the same number of records.
        Paths to single files will be uploaded as-is.
    session : `requests.Session`
        Connection to One Codex API.
    samples_resource : `onecodex.models.Samples`
        Wrapped potion-client object exposing `init_upload` and `confirm_upload` methods.
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
    A Sample object for the completed upload.
    """
    filename, file_size, file_format = _file_stats(files)

    # if filename cannot be represented as ascii, raise and suggest renaming
    try:
        # python2
        ascii_fname = unidecode(unicode(filename))
    except NameError:
        ascii_fname = unidecode(filename)

    if filename != ascii_fname:
        if coerce_ascii:
            # TODO: Consider warnings.warn here instead
            logging.warn(
                "Renaming {} to {}, must be ASCII\n".format(filename.encode("utf-8"), ascii_fname)
            )
            filename = ascii_fname
        else:
            raise OneCodexException("Filenames must be ascii. Try using --coerce-ascii")

    # disable progressbar while keeping context manager
    if not progressbar:
        progressbar = FakeProgressBar()

    # file_path is the path to the file on this disk. file_name is what we'll call the file in the
    # mainline database. file_size is the sum of both files in a pair, or the size of an unpaired
    # file. if paired, file_size is the uncompressed size. if unpaired, file_size is the actual
    # size on disk. unpaired files are uploaded as-is. paired files are decompressed, interleaved,
    # and uploaded as uncompressed data.
    with progressbar as bar:
        if isinstance(files, tuple):
            fobj = FASTXInterleave(files, file_size, file_format, bar)
        else:
            fobj = FilePassthru(files, file_size, bar)

        # must call init_upload in this loop in order to get a sample uuid we can call
        # cancel_upload on later if user hits ctrl+c
        fields = _call_init_upload(filename, file_size, metadata, tags, project, samples_resource)

        def cancel_atexit():
            bar.canceled = True
            bar.update(1)
            logging.info("Canceled upload for sample: {}".format(fields["sample_id"]))
            samples_resource.cancel_upload({"sample_id": fields["sample_id"]})

        atexit_register(cancel_atexit)

        # if the upload via init_upload fails, upload_sequence_fileobj will call
        # init_multipart_upload, which accepts metadata to be integrated into a newly-created
        # Sample model. if the s3 intermediate route is used, two Sample models will ultimately
        # exist on mainline: the failed fastx-proxy upload and the successful s3 intermediate.
        retry_fields = _make_retry_fields(filename, metadata, tags, project)

        try:
            sample_id = upload_sequence_fileobj(
                fobj, filename, fields, retry_fields, session, samples_resource
            )
            atexit_unregister(cancel_atexit)
            return sample_id
        except KeyboardInterrupt:
            cancel_atexit()
            atexit_unregister(cancel_atexit)
            raise


def _direct_upload(file_obj, file_name, fields, session, samples_resource):
    """Uploads a single file-like object via our validating proxy. Maintains compatibility with direct upload
    to a user's S3 bucket as well in case we disable our validating proxy.

    Parameters
    ----------
    file_obj : `FASTXInterleave`, `FilePassthru`, or a file-like object
        A wrapper around a pair of fastx files (`FASTXInterleave`) or a single fastx file. In the
        case of paired files, they will be interleaved and uploaded uncompressed. In the case of a
        single file, it will simply be passed through (`FilePassthru`) to One Codex, compressed
        or otherwise. If a file-like object is given, its mime-type will be sent as 'text/plain'.
    file_name : `string`
        The file_name you wish to associate this fastx file with at One Codex.
    fields : `dict`
        Additional data fields to include as JSON in the POST. Must include 'sample_id' and
        'upload_url' at a minimum.
    samples_resource : `onecodex.models.Samples`
        Wrapped potion-client object exposing `init_upload` and `confirm_upload` routes to mainline.

    Raises
    ------
    RetryableUploadException
        In cases where the proxy is temporarily down or we experience connectivity issues

    UploadException
        In other cases where the proxy determines the upload is invalid and should *not* be retried.
    """

    # need an OrderedDict to preserve field order for S3, required for Python 2.7
    multipart_fields = OrderedDict()

    for k, v in fields["additional_fields"].items():
        multipart_fields[str(k)] = str(v)

    # this attribute is only in FASTXInterleave and FilePassthru
    mime_type = getattr(file_obj, "mime_type", "text/plain")
    multipart_fields["file"] = (file_name, file_obj, mime_type)
    encoder = MultipartEncoder(multipart_fields)
    upload_request = None

    try:
        upload_request = session.post(
            fields["upload_url"],
            data=encoder,
            headers={"Content-Type": encoder.content_type},
            auth={},
        )
    except requests.exceptions.ConnectionError:
        pass

    # If we expect a status *always* try to check it,
    # waiting up to 15 minutes for buffering to complete
    if "status_url" in fields["additional_fields"]:
        now = time.time()
        while time.time() < (now + 60 * 15):
            try:
                resp = session.post(
                    fields["additional_fields"]["status_url"],
                    json={"sample_id": fields["sample_id"]},
                )
                resp.raise_for_status()
            except (ValueError, requests.exceptions.RequestException):
                raise RetryableUploadException(
                    "Unexpected failure of direct upload proxy. Retrying..."
                )

            if resp.json() and resp.json().get("complete", True) is False:
                logging.debug("Blocking on waiting for proxy to complete (in progress)...")
                time.sleep(5)
            else:
                break

        # Return is successfully processed
        if resp.json().get("code") in [200, 201]:
            file_obj.close()
            return
        elif resp.json().get("code") == 500:
            raise RetryableUploadException("Proxy failed. Retrying...")
        else:
            raise_api_error(resp, state="upload")

    # Direct to S3 case
    else:
        file_obj.close()
        if upload_request.status_code not in [200, 201]:
            raise RetryableUploadException("Unknown connectivity issue with proxy upload.")

        # Issue a callback -- this only happens in the direct-to-S3 case
        try:
            if not fields["additional_fields"].get("callback_url"):
                samples_resource.confirm_upload(
                    {"sample_id": fields["sample_id"], "upload_type": "standard"}
                )
        except requests.exceptions.HTTPError as e:
            raise_api_error(e.response, state="callback")
        except requests.exceptions.ConnectionError:
            raise_connectivity_error()


def upload_sequence_fileobj(file_obj, file_name, fields, retry_fields, session, samples_resource):
    """Uploads a single file-like object to the One Codex server via either fastx-proxy or directly
    to S3.

    Parameters
    ----------
    file_obj : `FASTXInterleave`, `FilePassthru`, or a file-like object
        A wrapper around a pair of fastx files (`FASTXInterleave`) or a single fastx file. In the
        case of paired files, they will be interleaved and uploaded uncompressed. In the case of a
        single file, it will simply be passed through (`FilePassthru`) to One Codex, compressed
        or otherwise. If a file-like object is given, its mime-type will be sent as 'text/plain'.
    file_name : `string`
        The file_name you wish to associate this fastx file with at One Codex.
    fields : `dict`
        Additional data fields to include as JSON in the POST. Must include 'sample_id' and
        'upload_url' at a minimum.
    retry_fields : `dict`
        Metadata sent to `init_multipart_upload` in the case that the upload via fastx-proxy fails.
    session : `requests.Session`
        Connection to One Codex API.
    samples_resource : `onecodex.models.Samples`
        Wrapped potion-client object exposing `init_upload` and `confirm_upload` routes to mainline.

    Raises
    ------
    UploadException
        In the case of a fatal exception during an upload.

    Returns
    -------
    `string` containing sample ID of newly uploaded file.
    """

    # First attempt to upload via our validating proxy
    try:
        _direct_upload(file_obj, file_name, fields, session, samples_resource)
        sample_id = fields["sample_id"]
    except RetryableUploadException:
        # upload failed--retry direct upload to S3 intermediate
        logging.error("{}: Connectivity issue, trying direct upload...".format(file_name))
        file_obj.seek(0)  # reset file_obj back to start

        try:
            retry_fields = samples_resource.init_multipart_upload(retry_fields)
        except requests.exceptions.HTTPError as e:
            raise_api_error(e.response, state="init")
        except requests.exceptions.ConnectionError:
            raise_connectivity_error(file_name)

        s3_upload = _s3_intermediate_upload(
            file_obj,
            file_name,
            retry_fields,
            session,
            samples_resource._client._root_url + retry_fields["callback_url"],  # full callback url
        )
        sample_id = s3_upload.get("sample_id", "<UUID not yet assigned>")

    logging.info("{}: finished as sample {}".format(file_name, sample_id))
    return sample_id


def upload_document(file_path, session, documents_resource, progressbar=None):
    """Uploads multiple document files to the One Codex server directly to S3 via an intermediate
    bucket.

    Parameters
    ----------
    file_path : `str`
        A path to a file on the system.
    session : `requests.Session`
        Connection to One Codex API.
    documents_resource : `onecodex.models.Documents`
        Wrapped potion-client object exposing `init_upload` and `confirm_upload` methods.
    progressbar : `click.progressbar`, optional
        If passed, display a progress bar using Click.

    Raises
    ------
    UploadException
        In the case of a fatal exception during an upload.

    Returns
    -------
    A `str` document ID for the newly uploaded file.
    """
    if not isinstance(file_path, six.string_types):
        raise ValueError(
            "Expected file_path to be a string, got {}".format(type(file_path).__name__)
        )

    file_name, file_size, _ = _file_stats(file_path, enforce_fastx=False)

    # disable progressbar while keeping context manager
    if not progressbar:
        progressbar = FakeProgressBar()

    with progressbar as bar:
        fobj = FilePassthru(file_path, file_size, bar)
        document_id = upload_document_fileobj(fobj, file_name, session, documents_resource)
        bar.finish()
        return document_id


def upload_document_fileobj(file_obj, file_name, session, documents_resource, log=None):
    """Uploads a single file-like object to the One Codex server directly to S3.

    Parameters
    ----------
    file_obj : `FilePassthru`, or a file-like object
        If a file-like object is given, its mime-type will be sent as 'text/plain'. Otherwise,
        `FilePassthru` will send a compressed type if the file is gzip'd or bzip'd.
    file_name : `string`
        The file_name you wish to associate this file with at One Codex.
    fields : `dict`
        Additional data fields to include as JSON in the POST.
    session : `requests.Session`
        Connection to One Codex API.
    documents_resource : `onecodex.models.Documents`
        Wrapped potion-client object exposing `init_upload` and `confirm_upload` routes to mainline.

    Notes
    -----
    In contrast to `upload_sample_fileobj`, this method will /only/ upload to an S3 intermediate
    bucket--not via our direct proxy or directly to a user's S3 bucket with a signed request.

    Raises
    ------
    UploadException
        In the case of a fatal exception during an upload.

    Returns
    -------
    `string` containing sample UUID of newly uploaded file.
    """
    try:
        fields = documents_resource.init_multipart_upload()
    except requests.exceptions.HTTPError as e:
        raise_api_error(e.response, state="init")
    except requests.exceptions.ConnectionError:
        raise_connectivity_error(file_name)

    s3_upload = _s3_intermediate_upload(
        file_obj,
        file_name,
        fields,
        session,
        documents_resource._client._root_url + fields["callback_url"],  # full callback url
    )

    document_id = s3_upload.get("document_id", "<UUID not yet assigned>")

    logging.info("{}: finished as document {}".format(file_name, document_id))
    return document_id


def _s3_intermediate_upload(file_obj, file_name, fields, session, callback_url):
    """Uploads a single file-like object to an intermediate S3 bucket which One Codex can pull from
    after receiving a callback.

    Parameters
    ----------
    file_obj : `FASTXInterleave`, `FilePassthru`, or a file-like object
        A wrapper around a pair of fastx files (`FASTXInterleave`) or a single fastx file. In the
        case of paired files, they will be interleaved and uploaded uncompressed. In the case of a
        single file, it will simply be passed through (`FilePassthru`) to One Codex, compressed
        or otherwise. If a file-like object is given, its mime-type will be sent as 'text/plain'.
    file_name : `string`
        The file_name you wish to associate this fastx file with at One Codex.
    fields : `dict`
        Additional data fields to include as JSON in the POST.
    callback_url : `string`
        API callback at One Codex which will trigger a pull from this S3 bucket.

    Raises
    ------
    UploadException
        In the case of a fatal exception during an upload. Note we rely on boto3 to handle its own retry logic.

    Returns
    -------
    `dict` : JSON results from internal confirm import callback URL
    """
    import boto3
    from boto3.s3.transfer import TransferConfig
    from boto3.exceptions import S3UploadFailedError

    # actually do the upload
    client = boto3.client(
        "s3",
        aws_access_key_id=fields["upload_aws_access_key_id"],
        aws_secret_access_key=fields["upload_aws_secret_access_key"],
    )

    # if boto uses threads, ctrl+c won't work
    config = TransferConfig(use_threads=False)

    # let boto3 update our progressbar rather than our FASTX wrappers, if applicable
    boto_kwargs = {}

    if hasattr(file_obj, "progressbar"):
        boto_kwargs["Callback"] = file_obj.progressbar.update
        file_obj.progressbar = None

    try:
        client.upload_fileobj(
            file_obj,
            fields["s3_bucket"],
            fields["file_id"],
            ExtraArgs={"ServerSideEncryption": "AES256"},
            Config=config,
            **boto_kwargs
        )
    except S3UploadFailedError:
        raise_connectivity_error(file_name)

    # issue a callback
    try:
        resp = session.post(
            callback_url,
            json={
                "s3_path": "s3://{}/{}".format(fields["s3_bucket"], fields["file_id"]),
                "filename": file_name,
                "import_as_document": fields.get("import_as_document", False),
            },
        )
    except requests.exceptions.ConnectionError:
        raise_connectivity_error(file_name)

    if resp.status_code != 200:
        raise_connectivity_error(file_name)

    try:
        return resp.json()
    except ValueError:
        return {}
