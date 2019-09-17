from __future__ import print_function, division

from collections import OrderedDict
import logging
import math
import requests
from requests_toolbelt import MultipartEncoder
import six
import time

from onecodex.exceptions import (
    OneCodexException,
    UploadException,
    RetryableUploadException,
    raise_connectivity_error,
    raise_api_error,
)
from onecodex.utils import atexit_register, atexit_unregister, snake_case
from onecodex.lib.files import FilePassthru, get_file_wrapper


log = logging.getLogger("onecodex")
DEFAULT_THREADS = 4


def _choose_boto3_chunksize(file_obj):
    """Return the appropriate chunksize for use in uploading the given file object.

    Choose the minimum chunk size for a boto3 direct-to-S3 upload that will result in less than
    10000 chunks (the maximum). This function will raise if there is no allowed chunk size big
    enough to accomodate the file.

    Parameters
    ----------
    file_obj : `FilePassthru`, or a file-like object
        In the case of a single file, it will simply be passed through (`FilePassthru`) to One Codex, compressed
        or otherwise.

    Returns
    -------
    `int`
        The minimum multipart chunk size in bytes.
    """
    file_obj_size = getattr(file_obj, "_fsize", None)

    if file_obj_size:
        allowed_chunk_sizes = [size * 1024 ** 2 for size in range(10, 110, 10)]

        for chunk_size in allowed_chunk_sizes:
            if math.ceil(file_obj_size / chunk_size) < 10000:
                break
        else:
            max_file_size = chunk_size * 10000

            raise OneCodexException(
                "File is too large to upload (size: {}, max: {})".format(
                    file_obj_size, max_file_size
                )
            )

        multipart_chunksize = chunk_size
    else:
        # default to 25 mb
        multipart_chunksize = 25 * 1024 ** 2

    return multipart_chunksize


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


def build_upload_dict(metadata, tags, project):
    """Build the metadata/tags/projects in a dict compatible with what the OneCodex backend expects."""
    upload_args = {}
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


def _call_init_upload(
    fobj, is_paired, metadata, tags, project, samples_resource, sample_id, external_sample_id
):
    """Call init_upload at the One Codex API and return data used to upload the file.

    Parameters
    ----------
    fobj : `PairedEndFiles` or `FilePassthru`
        A file object-like wrapper around the local file(s)
    is_paired: `bool`
    metadata : `dict`, optional
    tags : `list`, optional
    project : `string`, optional
        UUID of project to associate this sample with.
    samples_resource : `onecodex.models.Samples`
        Wrapped potion-client object exposing `init_upload` and `confirm_upload` routes to mainline.
    sample_id : `string`, optional
        If passed, will upload the file(s) to the sample with that id. Only works if the sample was pre-uploaded
    external_sample_id : `string`, optional
        If passed, will upload the file(s) to the sample with that metadata external id. Only works if the sample was pre-uploaded

    Returns
    -------
    `dict`
        Contains, at a minimum, 'upload_url' and 'sample_id'. Should also contain various additional
        data used to upload the file to fastx-proxy, a user's S3 bucket, or an intermediate bucket.
        If the input was some paired end files, it will also contain `paired_end_upload_url` and `paired_end_additional_fields` in
        order to upload the second file.
    """
    upload_args = {
        "upload_type": "standard",  # this is multipart form data
        "sample_id": sample_id,
        "external_sample_id": external_sample_id,
    }
    if is_paired:
        upload_args.update(
            {"filename": fobj.r1.filename, "paired_end_filename": fobj.r2.filename, "size": 1}
        )
    else:
        upload_args.update({"filename": fobj.filename, "size": fobj.size()})

    upload_args.update(build_upload_dict(metadata, tags, project))

    try:
        return samples_resource.init_upload(upload_args)
    except requests.exceptions.HTTPError as e:
        raise_api_error(e.response, state="init")
    except requests.exceptions.ConnectionError:
        raise_connectivity_error(upload_args["filename"])


def _make_retry_fields(
    fobj, is_paired, metadata, tags, project, sample_id=None, external_sample_id=None
):
    """Generate fields to send to init_multipart_upload.

    The fields returned by this function are used when a Sample upload via fastx-proxy fails.

    Parameters
    ----------
    fobj : `PairedEndFiles` or `FilePassthru`
        A file object-like wrapper around the local files
    is_paired: `bool`
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
    upload_args = {"sample_id": sample_id, "external_sample_id": external_sample_id}

    if is_paired:
        upload_args.update({"filename": fobj.r1.filename, "paired_end_filename": fobj.r2.filename})
    else:
        upload_args["filename"] = fobj.filename

    upload_args.update(build_upload_dict(metadata, tags, project))
    return upload_args


def preupload_sample(samples_resource, metadata=None, tags=None, project=None):
    """Make preupload request to the One Codex API and return the sample id.

    Parameters
    ----------
    metadata : `dict`, optional
    tags : `list`, optional
    project : `string`, optional
        UUID of project to associate this sample with.

    Returns
    -------
    `dict`
        Contains 'sample_id' field.
    """
    upload_args = build_upload_dict(metadata, tags, project)
    try:
        res = samples_resource.preupload(upload_args)
    except requests.exceptions.HTTPError as e:
        raise_api_error(e.response, state="init_preupload")

    return res["sample_id"]


def upload_sequence(
    file,
    session,
    samples_resource,
    metadata=None,
    tags=None,
    project=None,
    coerce_ascii=False,
    progressbar=None,
    sample_id=None,
    external_sample_id=None,
):
    """Upload a sequence file (or pair of files) to One Codex via our proxy or directly to S3.

    Parameters
    ----------
    file : `str` | `tuple(str, str)`
        A single file path or a tuple of paths for paired ends
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
    sample_id : `string`, optional
        If passed, will upload the file(s) to the sample with that id. Only works if the sample was pre-uploaded
    external_sample_id : `string`, optional
        If passed, will upload the file(s) to the sample with that metadata external id. Only works if the sample was pre-uploaded

    Returns
    -------
    A Sample object for the completed upload.
    """
    # disable progressbar while keeping context manager
    if not progressbar:
        progressbar = FakeProgressBar()

    # file_path is the path to the file on this disk. file_name is what we'll call the file in the
    # mainline database. file_size is the sum of both files in a pair, or the size of an unpaired
    # file. if paired, file_size is the uncompressed size. if unpaired, file_size is the actual
    # size on disk. unpaired files are uploaded as-is. paired files are decompressed, interleaved,
    # and uploaded as uncompressed data.
    with progressbar as bar:
        fobj = get_file_wrapper(file, coerce_ascii, bar)
        # So we don't have to check with isinstance which is going to be some mocks in tests
        is_paired = isinstance(file, tuple)

        # must call init_upload in this loop in order to get a sample uuid we can call
        # cancel_upload on later if user hits ctrl+c
        fields = _call_init_upload(
            fobj,
            is_paired,
            metadata,
            tags,
            project,
            samples_resource,
            sample_id,
            external_sample_id,
        )

        def cancel_atexit():
            bar.canceled = True
            bar.update(1)
            log.info("Canceled upload for sample: {}".format(fields["sample_id"]))

            try:
                samples_resource.cancel_upload({"sample_id": fields["sample_id"]})
            except requests.exceptions.HTTPError as e:
                # onecodex #298: it's possible to have this trigger after an upload has
                # already succeeded. try to catch that instead of blowing up
                if e.response and e.response.get("message") == "Upload already successful":
                    log.debug(
                        "Fail to cancel sample {}, upload already successful".format(
                            fields["sample_id"]
                        )
                    )
                else:
                    raise

        atexit_register(cancel_atexit)

        # if the upload via init_upload fails, upload_sequence_fileobj will call
        # init_multipart_upload, which accepts metadata to be integrated into a newly-created
        # Sample model. if the s3 intermediate route is used, two Sample models will ultimately
        # exist on mainline: the failed fastx-proxy upload and the successful s3 intermediate.
        # We also pass the `sample_id` and `external_sample_id`, which are typically None, to
        # support retries of pre-uploaded samples
        retry_fields = _make_retry_fields(
            fobj,
            is_paired,
            metadata,
            tags,
            project,
            sample_id=sample_id,
            external_sample_id=external_sample_id,
        )

        if is_paired:
            # 2 files to upload
            upload_sequence_fileobj(
                fobj.r1, fobj.r1.filename, fields, retry_fields, session, samples_resource
            )
            # TODO: check if we need to replace more and replace also in retry_fields
            fields["upload_url"] = fields["paired_end_upload_url"]
            fields["additional_fields"] = fields["paired_end_additional_fields"]
            sample_id = upload_sequence_fileobj(
                fobj.r2, fobj.r2.filename, fields, retry_fields, session, samples_resource
            )
        else:
            sample_id = upload_sequence_fileobj(
                fobj, fobj.filename, fields, retry_fields, session, samples_resource
            )

        atexit_unregister(cancel_atexit)
        return sample_id


def _direct_upload(file_obj, file_name, fields, session, samples_resource):
    """Upload a single file-like object via our validating proxy.

    Maintains compatibility with direct upload to a user's S3 bucket in case our validating proxy
    is disabled for this user.

    Parameters
    ----------
    file_obj : `FilePassthru`, or a file-like object
        In the case of a single file, it will simply be passed through (`FilePassthru`) to One Codex, compressed
        or otherwise. If a file-like object is given, its mime-type will be sent as 'text/plain'.
    file_name : `string`
        The file_name you wish to associate this fastx file with at One Codex.
    fields : `dict`
        Additional data fields to include as JSON in the POST. Must include 'sample_id' and
        'upload_url' at a minimum.
    session : `requests.Session`
        Use this session for direct uploads (via proxy or direct to a user's S3 bucket).
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

    # this attribute is only in FilePassthru
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
    # waiting up to 4 hours for buffering to complete (~30-50GB file gzipped)
    if "status_url" in fields["additional_fields"]:
        now = time.time()
        while time.time() < (now + 60 * 60 * 4):
            try:
                resp = session.post(
                    fields["additional_fields"]["status_url"],
                    json={"sample_id": fields["sample_id"]},
                )
                resp.raise_for_status()
            except (ValueError, requests.exceptions.RequestException) as e:
                log.debug("Retrying due to error: {}".format(e))
                raise RetryableUploadException(
                    "Unexpected failure of direct upload proxy. Retrying..."
                )

            if resp.json() and resp.json().get("complete", True) is False:
                log.debug("Blocking on waiting for proxy to complete (in progress)...")
                time.sleep(30)
            else:
                break

        # Return is successfully processed
        if resp.json().get("code") in [200, 201]:
            file_obj.close()
            return
        elif resp.json().get("code") == 500:
            log.debug("Retrying due to 500 from proxy...")
            raise RetryableUploadException("Unexpected issue with direct upload proxy. Retrying...")
        else:
            raise_api_error(resp, state="upload")

    # Direct to S3 case
    else:
        file_obj.close()
        if upload_request.status_code not in [200, 201]:
            raise UploadException("Unknown connectivity issue with direct upload.")

        # Issue a callback -- this only happens in the direct-to-S3 case
        try:
            if not fields["additional_fields"].get("callback_url"):
                samples_resource.confirm_upload(
                    {"sample_id": fields["sample_id"], "upload_type": "standard"}
                )
        except requests.exceptions.HTTPError as e:
            raise_api_error(e.response, state="callback")
        except requests.exceptions.ConnectionError:
            raise_connectivity_error(file_name)


def upload_sequence_fileobj(file_obj, file_name, fields, retry_fields, session, samples_resource):
    """Upload a single file-like object to One Codex via either fastx-proxy or directly to S3.

    Parameters
    ----------
    file_obj : `FilePassthru`, or a file-like object
        In the case of a single file, it will simply be passed through (`FilePassthru`) to One Codex, compressed
        or otherwise. If a file-like object is given, its mime-type will be sent as 'text/plain'.
    file_name : `string`
        The file_name you wish to associate this fastx file with at One Codex.
    fields : `dict`
        Additional data fields to include as JSON in the POST. Must include 'sample_id' and
        'upload_url' at a minimum.
    retry_fields : `dict`
        Metadata sent to `init_multipart_upload` in the case that the upload via fastx-proxy fails.
    session : `requests.Session`
        Use this session for direct uploads (via proxy or direct to a user's S3 bucket).
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
        sample_id = fields["sample_id"]

        # Are we being directed to skip the proxy? If so, only do it if files ares <5GB since that's the limit for
        # direct uploads to S3
        if (
            "AWSAccessKeyId" in fields["additional_fields"]
            and getattr(file_obj, "_fsize", 0) > 5 * 1024 ** 3
        ):
            raise RetryableUploadException

        # Big files are going to skip the proxy even if the backend told us the opposite
        # 100GB is considered big enough to defer the validation
        # In some cases, file_obj might be a BytesIO object instead of one of our file object so we
        # filter them out by checking for a `write` attribute
        if not hasattr(file_obj, "write") and file_obj.size() > 100 * 1024 ** 3:
            raise RetryableUploadException

        _direct_upload(file_obj, file_name, fields, session, samples_resource)
    except RetryableUploadException:
        # upload failed -- retry direct upload to S3 intermediate; first try to cancel pending upload
        try:
            samples_resource.cancel_upload({"sample_id": sample_id})
        except Exception as e:
            log.debug("Failed to cancel upload: {}".format(e))
        log.error("{}: Connectivity issue, trying upload via intermediary...".format(file_name))
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
            samples_resource._client.session,
            samples_resource._client._root_url + retry_fields["callback_url"],  # full callback url
        )
        sample_id = s3_upload.get("sample_id", "<UUID not yet assigned>")

    log.info("{}: finished as sample {}".format(file_name, sample_id))
    return sample_id


def upload_document(file_path, documents_resource, progressbar=None):
    """Upload multiple document files to One Codex directly to S3 via an intermediate bucket.

    Parameters
    ----------
    file_path : `str`
        A path to a file on the system.
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

    # disable progressbar while keeping context manager
    if not progressbar:
        progressbar = FakeProgressBar()

    with progressbar as bar:
        fobj = FilePassthru(file_path, bar)
        document_id = upload_document_fileobj(fobj, fobj.filename, documents_resource)
        bar.finish()
        return document_id


def upload_document_fileobj(file_obj, file_name, documents_resource):
    """Upload a single file-like object to One Codex directly to S3 via an intermediate bucket.

    Parameters
    ----------
    file_obj : `FilePassthru`, or a file-like object
        If a file-like object is given, its mime-type will be sent as 'text/plain'. Otherwise,
        `FilePassthru` will send a compressed type if the file is gzip'd or bzip'd.
    file_name : `string`
        The file_name you wish to associate this file with at One Codex.
    fields : `dict`
        Additional data fields to include as JSON in the POST.
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
        documents_resource._client.session,
        documents_resource._client._root_url + fields["callback_url"],  # full callback url
    )

    document_id = s3_upload.get("document_id", "<UUID not yet assigned>")

    log.info("{}: finished as document {}".format(file_name, document_id))
    return document_id


def _s3_intermediate_upload(file_obj, file_name, fields, session, callback_url):
    """Upload a single file-like object to an intermediate S3 bucket.

    One Codex will pull the file from S3 after receiving a callback.

    Parameters
    ----------
    file_obj : `FilePassthru`, or a file-like object
        In the case of a single file, it will simply be passed through (`FilePassthru`) to One Codex, compressed
        or otherwise. If a file-like object is given, its mime-type will be sent as 'text/plain'.
    file_name : `string`
        The file_name you wish to associate this fastx file with at One Codex.
    fields : `dict`
        Additional data fields to include as JSON in the POST.
    session : `requests.Session`
        Authenticated connection to One Codex API used to POST callback.
    callback_url : `string`
        API callback at One Codex which will trigger a pull from this S3 bucket.

    Raises
    ------
    UploadException
        In the case of a fatal exception during an upload. Note we rely on boto3 to handle its own
        retry logic.

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

    multipart_chunksize = _choose_boto3_chunksize(file_obj)

    # if boto uses threads, ctrl+c won't work
    config = TransferConfig(use_threads=False, multipart_chunksize=multipart_chunksize)

    # let boto3 update our progressbar rather than our FASTX wrappers, if applicable
    boto_kwargs = {}

    if hasattr(file_obj, "progressbar"):
        boto_kwargs["Callback"] = file_obj.progressbar.update
        file_obj._progressbar = file_obj.progressbar
        file_obj.progressbar = None

    for attempt in range(1, 4):
        try:
            client.upload_fileobj(
                file_obj,
                fields["s3_bucket"],
                fields["file_id"],
                ExtraArgs={"ServerSideEncryption": "AES256"},
                Config=config,
                **boto_kwargs
            )
            break
        except S3UploadFailedError as e:
            logging.debug("Caught S3UploadFailedError on attempt {}/3: {}".format(attempt, str(e)))
            logging.error(
                "{}: Connectivity issue, retrying upload via intermediary ({}/3)...".format(
                    file_name, attempt
                )
            )

            # rewind the progressbar if possible, then remove so boto3 can update the bar directly
            if hasattr(file_obj, "_progressbar"):
                file_obj.progressbar = file_obj._progressbar
                file_obj.seek(0)
                file_obj.progressbar = None
            else:
                file_obj.seek(0)
    else:
        logging.debug("{}: exhausted all retries via intermediary")
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
