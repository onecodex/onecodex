from __future__ import print_function, division

import atexit
import copy
import logging
import math
import requests
import six

from onecodex.exceptions import OneCodexException, raise_connectivity_error, raise_api_error
from onecodex.utils import snake_case, FakeProgressBar
from onecodex.lib.files import FilePassthru, get_file_wrapper

from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry


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
        allowed_chunk_sizes = [size * 1024**2 for size in range(10, 110, 10)]

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
        multipart_chunksize = 25 * 1024**2

    return multipart_chunksize


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


def _get_init_multipart_upload_payload(
    fobj, is_paired, metadata, tags, project, sample_id=None, external_sample_id=None
):
    """Generate fields to send to init_multipart_upload.

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
    samples_resource,
    metadata=None,
    tags=None,
    project=None,
    coerce_ascii=False,
    progressbar=None,
    sample_id=None,
    external_sample_id=None,
):
    """Upload a sequence file (or pair of files) to One Codex directly to S3.

    Parameters
    ----------
    file : `str` | `tuple(str, str)`
        A single file path or a tuple of paths for paired ends
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

    with progressbar as bar:
        fobj = get_file_wrapper(file, coerce_ascii, bar)
        # So we don't have to check with isinstance which is going to be some mocks in tests
        is_paired = isinstance(file, tuple)

        # must call init_multipart_upload in this loop in order to get a sample uuid we can call
        # cancel_upload on later if user hits ctrl+c

        # upload_sequence_fileobj will call init_multipart_upload, which accepts metadata to be integrated
        # into a newly-created Sample model.
        # We also pass the `sample_id` and `external_sample_id`, which are typically None, to
        # support retries of pre-uploaded samples
        payload = _get_init_multipart_upload_payload(
            fobj,
            is_paired,
            metadata,
            tags,
            project,
            sample_id=sample_id,
            external_sample_id=external_sample_id,
        )

        try:
            fields = samples_resource.init_multipart_upload(payload)
        except requests.exceptions.HTTPError as e:
            raise_api_error(e.response, state="init")

        def cancel_atexit():
            bar.canceled = True
            bar.update(1)

            if is_paired:
                filename = "{} and {}".format(fobj.r1.filename, fobj.r2.filename)
            else:
                filename = fobj.filename

            log.info("Canceled upload for {} as sample {}".format(filename, fields["sample_id"]))

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

        atexit.register(cancel_atexit)

        if is_paired:
            # 2 files to upload
            # The backend will check for the r1 file in the callback so we upload r2 first
            fields_pe = copy.deepcopy(fields)
            fields_pe["file_id"] = fields_pe["paired_end_file_id"]
            _upload_sequence_fileobj(
                fobj.r2, fobj.r2.filename, fields_pe, samples_resource, callback=False
            )
            sample_id = _upload_sequence_fileobj(
                fobj.r1, fobj.r1.filename, fields, samples_resource
            )
        else:
            sample_id = _upload_sequence_fileobj(fobj, fobj.filename, fields, samples_resource)

        atexit.unregister(cancel_atexit)
        return sample_id


def _upload_sequence_fileobj(file_obj, file_name, fields, samples_resource, callback=True):
    """Upload a single file-like object to One Codex to S3.

    Parameters
    ----------
    file_obj : `FilePassthru`, or a file-like object
        In the case of a single file, it will simply be passed through (`FilePassthru`) to One Codex, compressed
        or otherwise. If a file-like object is given, its mime-type will be sent as 'text/plain'.
    file_name : `string`
        The file_name you wish to associate this fastx file with at One Codex.
    fields : `dict`
        The fields boto will need to have to upload to S3
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
    s3_upload = _s3_intermediate_upload(
        file_obj,
        file_name,
        fields,
        samples_resource._client.session,
        samples_resource._client._root_url + fields["callback_url"]
        if callback
        else None,  # full callback url
    )
    sample_id = s3_upload.get("sample_id")

    msg = "{}: upload finished".format(file_name)
    if sample_id is not None:
        msg += " as sample {}".format(sample_id)

    log.info(msg)
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
        document_id = _upload_document_fileobj(fobj, fobj.filename, documents_resource)
        bar.finish()
        return document_id


def _upload_document_fileobj(file_obj, file_name, documents_resource):
    """Upload a single file-like object to One Codex directly to S3 via an intermediate bucket.

    Parameters
    ----------
    file_obj : `FilePassthru`, or a file-like object
        If a file-like object is given, its mime-type will be sent as 'text/plain'. Otherwise,
        `FilePassthru` will send a compressed type if the file is gzip'd or bzip'd.
    file_name : `string`
        The file_name you wish to associate this file with at One Codex.
    documents_resource : `onecodex.models.Documents`
        Wrapped potion-client object exposing `init_upload` and `confirm_upload` routes to mainline.

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

    msg = "{}: finished".format(file_name)
    document_id = s3_upload.get("document_id")
    if document_id is not None:
        msg += " as document {}".format(document_id)

    log.info(msg)
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

    boto3_session = boto3.session.Session()
    # actually do the upload
    client = boto3_session.client(
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

    # In paired uploads, we only want to call the callback url once both files are uploaded
    if not callback_url:
        return {}

    # issue a callback
    try:
        # retry on 502, 503, 429, with a backoff timing of 4s, 8s, and 16s, False retries on all HTTP methods
        retry_strategy = Retry(
            total=3, backoff_factor=4, method_whitelist=False, status_forcelist=[502, 503, 429]
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount(callback_url, adapter)

        resp = session.post(
            callback_url,
            json={
                "s3_path": "s3://{}/{}".format(fields["s3_bucket"], fields["file_id"]),
                "filename": file_name,  #
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
