import json
import logging
from multiprocessing import Lock, Value
import os
import sys
from threading import BoundedSemaphore, Thread

import requests
from requests_toolbelt import MultipartEncoder, MultipartEncoderMonitor

from onecodex.utils import check_for_allowed_file, warn_if_insecure_platform


# Use multipart upload for anything over ~5GB
# (a little bit smaller bc 1000^3, not 1024^3)
MULTIPART_SIZE = 5 * 1000 * 1000 * 1000
DEFAULT_UPLOAD_THREADS = 4
CHUNK_SIZE = 8192

log = logging.getLogger(__name__)


def old_upload(files, session, samples_resource, server_url, threads=None):
    """
    This is the entry point for the upload flow. It will determine
        what approach to take with uploading and pass to other functions

    -param list files: The list of file (paths) to upload
    -param int threads: Number of upload threads to use (def=4)
    """
    if threads is None:
        threads = DEFAULT_UPLOAD_THREADS

    # check insecure platform, disable warnigns
    if warn_if_insecure_platform():
        logging.captureWarnings(True)

    file_sizes = [os.path.getsize(f) for f in files]
    if min(file_sizes) < 35:
        print("Cannot upload empty files. Please check that all files "
              "contain sequence data and try again.")
        raise SystemExit

    max_file_size = max(file_sizes)
    if max_file_size > MULTIPART_SIZE:
        for ix, f in enumerate(files):
            if file_sizes[ix] > MULTIPART_SIZE:
                _upload_multipart(f, session, samples_resource, server_url)
            else:
                _upload_direct([f], session, samples_resource, server_url, threads)
    else:
        _upload_direct(files, session, samples_resource, server_url, threads)


def _upload_multipart(filename, session, samples_resource, server_url):
    """
    This is the upload function for files over 5GB. It uploads them serially
        to s3 using awscli. It will exit if awscli is not installed

    -param str file: The path to the file to upload
    """
    check_for_allowed_file(filename)
    multipart_req = samples_resource.read_init_multipart_upload()

    s3_bucket = multipart_req["s3_bucket"]
    callback_chunk = multipart_req['callback_url']
    callback_url = server_url.rstrip("/") + callback_chunk
    file_id = multipart_req["file_id"]
    aws_access_key_id = multipart_req["upload_aws_access_key_id"]
    aws_secret_access_key = multipart_req["upload_aws_secret_access_key"]  # noqa

    # Upload to s3 using boto
    try:
        import awscli  # noqa
        import subprocess
    except ImportError:
        print("You must install the awscli package for files >5GB in size. "
              "On most systems, it can be installed with `pip install awscli`.")
        raise SystemExit

    s3_path = "s3://" + s3_bucket + "/" + file_id
    print("Starting large (>5GB) file upload. "
          "Please be patient while the file transfers...")
    try:
        # We want to only get output from onecodex
        p = subprocess.Popen("AWS_ACCESS_KEY_ID=%s AWS_SECRET_ACCESS_KEY=%s aws s3 cp %s %s --sse" %
                             (aws_access_key_id, aws_secret_access_key, filename, s3_path),
                             shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        print("\n"
              "    ###########################################################\n"
              "    ###           Uploading large multipart file            ###\n"
              "    ###                Upload output below...               ###\n"
              "    ###########################################################\n")
        while p.poll() is None:
            char = p.stdout.read(1)
            sys.stdout.write(str(char))
            sys.stdout.flush()

    except KeyboardInterrupt:
        log.info("Upload successfully cancelled. Quitting.")
        p.sigterm()
        raise SystemExit

    if p.returncode != 0:
        print('An error occured uploading %s using the aws-cli.' % filename)
        raise SystemExit

    callback_request = session.post(callback_url,
                                    headers={"Content-Type": "application/json"},
                                    data=json.dumps({"s3_path": s3_path,
                                                     "filename": os.path.basename(filename)}))

    if callback_request.status_code != 200:
        print("Upload of %s failed. Please contact help@onecodex.com "
              "if you experience further issues." % filename)
        sys.exit(1)
    print("Successfully uploaded: %s\n" % filename)
    print("    ###########################################################\n"
          "    ### Please note: Large file uploads may take several    ###\n"
          "    ### minutes to appear on the One Codex website. If a    ###\n"
          "    ### file does not appear after a longer period of time, ###\n"
          "    ### however, please contact us at help@onecodex.com.    ###\n"
          "    ###########################################################\n")


def _upload_direct(files, session, samples_resource, server_url, threads):
    """
    This is the upload method for files < 5GB. They are sent directly
        with requests and using the threads

    -param list files: The list of filepaths to upload
    -param int threads: Number of threads to use
    """

    if threads > 1:
        semaphore = BoundedSemaphore(threads)
    if threads != DEFAULT_UPLOAD_THREADS:
        print("Uploading with up to %d thread(s)." % threads)

    # Set up upload threads
    upload_threads = []
    upload_progress_bytes = Value('L', 0)
    upload_progress_lock = Lock()
    total_bytes = sum([os.path.getsize(f) for f in files])
    total_files = Value('i', len(files))

    for f in files:
        if threads > 1 and len(files) > 1:  # parallel uploads
            # Multi-threaded uploads
            t = Thread(target=_upload_helper,
                       args=(f, session, samples_resource, server_url,
                             upload_progress_bytes, upload_progress_lock,
                             total_bytes, total_files, semaphore))
            upload_threads.append(t)
            t.start()
        else:  # serial uploads
            _upload_helper(f, session, samples_resource, server_url,
                           upload_progress_bytes, upload_progress_lock,
                           total_bytes, total_files)

        if threads > 1:
            for ut in upload_threads:
                ut.join()


def _upload_helper(filename, session, samples_resource, server_url,
                   upload_progress_bytes, upload_progress_lock,
                   total_bytes, total_files, semaphore=None):
    """
    This is the tread worker function for direct uploads. It makes several
        calls to the app server to sign the upload and record success.
        It also passes byte amount info to the callback for prograss bar

    -param str filename: The filepath to be uploaded in this thread
    -param Resource samples_resource: The Sample API resource
    -param str server_url: The server's base URL
    -param int upload_progress_bytes: Bytes uploaded so far
    -param Lock upload_progress_lock: The thread lock
    -param int total_bytes: Total bytes left to upload
    -param int total_files: Total number of files to upload
    -param BoundedSemaphore semaphore: Count of threads in existence

    """
    stripped_filename = os.path.basename(filename)
    try:
        upload_info = samples_resource.init_upload({
            'filename': stripped_filename,
            'size': os.path.getsize(filename),
            'upload_type': 'standard'  # This is multipart form data
        })
    except requests.exceptions.HTTPError:
        print('The attempt to initiate your upload failed. Please make '
              'sure you are logged in (`onecodex login`) and try again. '
              'If you continue to experience problems, contact us at '
              'help@onecodex.com for assistance.')
        raise SystemExit
    upload_url = upload_info['upload_url']

    # First get the signing form data
    if semaphore is not None:
        semaphore.acquire()

    # Coerce to str or MultipartEncoder fails
    # Need a list to preserve order for S3
    fields = []
    for k, v in upload_info['additional_fields'].items():
        fields.append((str(k), str(v)))

    fields.append(("file", (stripped_filename, open(filename, mode='rb'), "text/plain")))
    e = MultipartEncoder(fields)
    m = MultipartEncoderMonitor(e, lambda x: _upload_callback(x, upload_progress_bytes,
                                                              upload_progress_lock,
                                                              total_bytes=(total_bytes + 8192),
                                                              n_files=total_files))

    max_retries = 3
    n_retries = 0
    while n_retries < max_retries:
        try:
            upload_request = session.post(upload_url, data=m,
                                          headers={"Content-Type": m.content_type},
                                          auth={})
            if upload_request.status_code != 201:
                print("Upload failed. Please contact help@onecodex.com for assistance.")
                raise SystemExit
            break
        except requests.exceptions.ConnectionError:
            n_retries += 1
            if n_retries == max_retries:
                print("The command line client is experiencing connectivity issues and "
                      "cannot complete the upload of %s at this time. Please try again "
                      "later. If the problem persists, contact us at help@onecodex.com "
                      "for assistance." % stripped_filename)
                raise SystemExit

    # Finally, issue a callback
    try:
        samples_resource.confirm_upload({
            'sample_id': upload_info['sample_id'],
            'upload_type': 'standard'
        })
        success_msg = ("Successfully uploaded: %s. Sample ID is: %s." %
                       (filename, upload_info['sample_id']))
        if upload_progress_bytes.value == -1:  # == -1 upon completion
            print(success_msg)
        else:
            sys.stderr.write("\r")
            sys.stderr.flush()
            print(success_msg)
        with upload_progress_lock:
            total_files.value -= 1
    except requests.exceptions.HTTPError:
        print("Failed to upload: %s" % filename)
        raise SystemExit

    if semaphore is not None:
        semaphore.release()


def _upload_callback(monitor, upload_progress_bytes, lock, total_bytes, n_files):
    """
    This is the callback/monitor function for the upload threads. It uses
        the byte information to make a progress bar and what not.
        -param int monitor
        -param int upload_progress_bytes
        -param Lock lock
        -param int total_bytes
        -param int n_files
    """
    if upload_progress_bytes.value == -1:
        return
    with lock:
        upload_progress_bytes.value += CHUNK_SIZE  # Chunk size
    if upload_progress_bytes.value == 0:
        progress = 0.0
    else:
        progress = upload_progress_bytes.value / float(total_bytes)
    bar_length = 20
    if progress < 0:
        progress = 0
        status = "Halt...                       \r\n"
    elif progress >= 1:
        progress = 1
        status = "Done.                         \r\n"
        with lock:
            upload_progress_bytes.value = -1
    elif progress <= 0.9:
        status = ''
    else:
        status = 'Completing...'
    block = int(round(bar_length * progress))
    text = "\rUploading: [{0}] {1:.2f}% {2}".format(
        "#" * block + "-" * (bar_length - block),
        progress * 100, status)

    sys.stderr.write(text)
    sys.stderr.flush()
