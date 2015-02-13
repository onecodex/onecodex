"""
Functions to implement the v0 One Codex API calls.
"""
from __future__ import print_function
import json
import os
import requests
from requests_toolbelt import MultipartEncoder, MultipartEncoderMonitor
import sys
from multiprocessing import Lock, Value
from threading import BoundedSemaphore, Thread
import urlparse
from onecodex import version
from onecodex.helpers import stderr


# Config
if os.environ.get("ONE_CODEX_API_BASE") is not None:
    BASE_API = os.environ.get("ONE_CODEX_API_BASE")
    print("ALL REQUESTS GOING THROUGH: %s" % BASE_API)
else:
    BASE_API = "https://beta.onecodex.com/api/v0/"

BASE_URL = urlparse.urlparse(BASE_API)
BASE_URL = BASE_URL._replace(path='/').geturl()
DEFAULT_THREADS = 4
CHUNK_SIZE = 8192

# Use multipart upload for anything over ~5GB
# (a little bit smaller bc 1000^3, not 1024^3)
MULTIPART_SIZE = 5 * 1000 * 1000 * 1000

BAD_AUTH_MSG = ("\nYour login credentials appear be bad. Try logging out:"
                "\n    onecodex logout"
                "\n"
                "\nAnd then logging back in:"
                "\n    onecodex login"
                "\n")

BAD_API_KEY_MSG = ("\nThe --api-key you entered appears to be "
                   "invalid. Please double check the key and try again.\n")


# Helpers
def pprint(j, args):
    if args.pprint:
        print(json.dumps(j, sort_keys=True,
                         indent=4, separators=(',', ': ')))
    else:
        print(j)


def download_file_helper(url, input_path, auth=None):
    r = requests.get(url, stream=True, auth=auth)
    if r.status_code != 200:
        stderr("Failed to download file: %s" % r.json()["message"])
    original_filename = urlparse.urlparse(r.url).path.split("/")[-1]
    if os.path.isdir(input_path):
        local_full_path = os.path.join(input_path, original_filename)
    else:
        local_full_path = input_path
    with open(local_full_path, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024):
            if chunk:  # filter out keep-alive new chunks
                f.write(chunk)
                f.flush()
    print("Successfully downloaded %s to %s" % (original_filename, local_full_path))


# Version checking function
def get_update_message():
    r = requests.post(BASE_API + "check_for_cli_update",
                      data={"version": version.VERSION,
                            "api_version": version.API_VERSION})
    if r.status_code == 200:
        j = r.json()
        if j.get("message"):
            stderr(j["message"])


def upload_callback(monitor, upload_progress_bytes, lock, total_bytes, n_files):
    if upload_progress_bytes.value == -1:
        return
    with lock:
        upload_progress_bytes.value += CHUNK_SIZE  # Chunk size
    if upload_progress_bytes.value == 0:
        progress = 0.0
    else:
        progress = upload_progress_bytes.value / float(total_bytes)
    barLength = 20  # Modify this to change the length of the progress bar
    if progress < 0:
        progress = 0
        status = "Halt...                       \r\n"  # Needs to be longer than files string
    elif progress >= 1:
        progress = 1
        status = "Done.                         \r\n"
        with lock:
            upload_progress_bytes.value = -1
    else:
        status = ("(%d files remaining)" % n_files.value
                  if n_files.value > 1
                  else "(1 file remaining)")
    block = int(round(barLength * progress))
    text = "\rUploading: [{0}] {1:.2f}% {2}".format("#" * block + "-" * (barLength - block),
                                                    progress * 100, status)
    sys.stderr.write(text)
    sys.stderr.flush()


# Upload functions
def upload(args):
    """
    Note that this doesn't actually use the default API route -- it instead
    posts directly to S3.
    """
    file_sizes = [os.path.getsize(f) for f in args.file]
    max_file_size = max(file_sizes)
    if max_file_size > MULTIPART_SIZE:
        for ix, f in enumerate(args.file):
            if file_sizes[ix] > MULTIPART_SIZE:
                upload_multipart(args, f)
            else:
                upload_direct(args, [f])
    else:
        upload_direct(args, args.file)


def upload_multipart(args, f):
    """
    Note, for large files we upload them one at a time
    using a special API.
    """
    creds = (args.credentials['api_key'], '')
    r0 = requests.get(BASE_API + "init_multipart_upload", auth=creds)
    if r0.status_code != 200:
        stderr("Failed to initiate large multipart upload (>5GB).")
        sys.exit(1)

    s3_bucket = r0.json()["s3_bucket"]
    callback_url = BASE_URL.rstrip("/") + r0.json()['callback_url']
    file_id = r0.json()["file_id"]

    # Upload to s3 using boto
    try:
        import awscli  # noqa
        import subprocess
    except ImportError:
        stderr("You must install the awscli package for files >5GB in size. "
               "On most systems, it can be installed with `pip install awscli`.")
        sys.exit(1)

    s3_path = "s3://" + s3_bucket + "/" + file_id
    print("Starting large (>5GB) file upload. Please be patient while the file transfers...")
    try:
        p = subprocess.Popen("aws s3 cp %s %s" % (f, s3_path),
                             stderr=subprocess.STDOUT, shell=True)
        p.wait()
    except KeyboardInterrupt:
        print("Upload successfully cancelled. Quitting.")
        p.sigterm()
        sys.exit(1)

    if p.returncode != 0:
        stderr("Failed to upload %s" % f)
        sys.exit(1)

    r1 = requests.post(callback_url, auth=creds,
                       headers={"Content-Type": "application/json"},
                       data=json.dumps({"s3_path": s3_path,
                                        "filename": os.path.basename(f)}))
    if r1.status_code != 200:
        stderr("Upload of %s failed. Please contact help@onecodex.com "
               "if you experience further issues." % f)
        sys.exit(1)
    print("Successfully uploaded: %s\n" % f)
    print("    ###########################################################\n"
          "    ### Please note: Large file uploads may take several    ###\n"
          "    ### minutes to appear on the One Codex website. If a    ###\n"
          "    ### file does not appear after a longer period of time, ###\n"
          "    ### however, please contact us at help@onecodex.com.    ###\n"
          "    ###########################################################\n")


def upload_direct(args, files):
    """
    Directly POST to S3. This does not use s3cmd or multipart uploads.
    """
    creds = (args.credentials['api_key'], '')

    if args.threads:
        semaphore = BoundedSemaphore(args.max_threads)
        if args.max_threads != DEFAULT_THREADS:
            print("Uploading with up to %d threads." % args.max_threads)

    # Get the initially needed routes
    r0 = requests.get(BASE_API + 'presign_upload', auth=creds)
    if r0.status_code == 401:
        if not args.api_key:
            stderr(BAD_AUTH_MSG)
        else:
            stderr(BAD_API_KEY_MSG)
        sys.exit(1)
    elif r0.status_code != 200:
        stderr("Failed to get upload signing credentials")
        sys.exit(1)

    j0 = r0.json()
    s3_url = j0['url']
    signing_url = BASE_URL.rstrip("/") + j0['signing_url']
    callback_url = BASE_URL.rstrip("/") + j0['callback_url']

    upload_threads = []
    upload_progress_bytes = Value('L', 0)
    upload_progress_lock = Lock()
    total_bytes = sum([os.path.getsize(f) for f in files])
    total_files = Value('i', len(files))
    for f in files:
        if args.threads and len(files) > 1:  # parallel uploads
            # Multi-threaded uploads
            t = Thread(target=upload_helper,
                       args=(f, s3_url, signing_url, callback_url,
                             creds, upload_progress_bytes, upload_progress_lock,
                             total_bytes, total_files, semaphore))
            upload_threads.append(t)
            t.start()
        else:  # serial uploads
            upload_helper(f, s3_url, signing_url, callback_url, creds,
                          upload_progress_bytes, upload_progress_lock,
                          total_bytes, total_files)

    if args.threads:
        for ut in upload_threads:
            ut.join()


def upload_helper(f, s3_url, signing_url, callback_url, creds,
                  upload_progress_bytes, upload_progress_lock,
                  total_bytes, total_files,
                  semaphore=None):
    # First get the signing form data
    if semaphore is not None:
        semaphore.acquire()

    stripped_filename = os.path.basename(f)
    r1 = requests.post(signing_url, data={"filename": stripped_filename, "via_api": "true"},
                       auth=creds)
    if r1.status_code != 200:
        try:
            stderr("Failed upload: %s" % r1.json()["msg"])
        except:
            stderr("Upload failed. Please contact help@onecodex.com for "
                   "assistance if you continue to experience problems.")
        sys.exit(1)
    file_uuid = r1.json()['key'].split("/")[-2][5:]

    # Coerce to str or MultipartEncoder fails
    # Need a list to preserve order for S3
    fields = []
    for k, v in r1.json().items():
        fields.append((str(k), str(v)))

    fields.append(("file", (stripped_filename, open(f, mode='rb'), "text/plain")))
    e = MultipartEncoder(fields)
    m = MultipartEncoderMonitor(e, lambda x: upload_callback(x, upload_progress_bytes,
                                                             upload_progress_lock,
                                                             total_bytes=(total_bytes + 8192),
                                                             n_files=total_files))
    max_retries = 3
    n_retries = 0
    while n_retries < max_retries:
        try:
            r2 = requests.post(s3_url, data=m, headers={"Content-Type": m.content_type})
            if r2.status_code != 201:
                stderr("Upload failed. Please contact help@onecodex.com for assistance.")
                sys.exit(1)
            break
        except requests.exceptions.ConnectionError:
            n_retries += 1
            if n_retries == max_retries:
                stderr("The command line client is experiencing connectivity issues and "
                       "cannot complete the upload of %s at this time. Please try again "
                       "later. If the problem persists, contact us at help@onecodex.com "
                       "for assistance." % stripped_filename)
                sys.exit(1)

    # Finally, issue a callback
    r3 = requests.post(callback_url, auth=creds, data={
        "location": r2.headers['location'],
        "size": os.path.getsize(f)
    })
    if r3.status_code == 200:
        success_msg = "Successfully uploaded: %s. File ID is: %s." % (f, file_uuid)
        if upload_progress_bytes.value == -1:  # == -1 upon completion
            print(success_msg)
        else:
            sys.stderr.write("\r")
            sys.stderr.flush()
            print(success_msg)
        with upload_progress_lock:
            total_files.value -= 1
    else:
        print("Failed to upload: %s" % f)
        sys.exit(1)

    if semaphore is not None:
        semaphore.release()


# Helper for /route/UUID pattern
def api_helper(args, route, supplement=""):
    creds = (args.credentials['api_key'], '')
    if not getattr(args, route):
        r = requests.get(BASE_API + route + supplement,
                         auth=creds)
        abort_helper(r, args)
        j = r.json()
        pprint(j, args)
    else:
        for uuid in getattr(args, route):
            r = requests.get(BASE_API + route + "/" + uuid + supplement,
                             auth=creds)
            abort_helper(r, args)
            j = r.json()
            pprint(j, args)


def abort_helper(r, args):
    if r.status_code == 401:
        if not args.api_key:
            stderr(BAD_AUTH_MSG)
        else:
            stderr(BAD_API_KEY_MSG)
        sys.exit(1)


def samples(args):
    api_helper(args, route="samples")


def analyses(args):
    if not args.raw and not args.table:
        api_helper(args, route="analyses")
    elif args.raw and args.table:
        stderr("Can only request raw or table data at the same time.")
        sys.exit(1)
    elif args.raw and not args.table:
        if len(args.analyses) == 0:
            stderr("No analysis specified. Please note the first argument "
                   "following --raw is an optional path for storing the raw "
                   "download.\nIf you do not want to specify this, append --raw "
                   "to the end of your command, e.g., `onecodex analyses <id> --raw`.")
            sys.exit(1)
        elif len(args.analyses) != 1:
            stderr("Can only request raw data on one Analysis at a time.")
            sys.exit(1)
        download_file_helper(BASE_API + "analyses/" + args.analyses[0] + "/raw",
                             input_path=args.raw,
                             auth=(args.credentials['api_key'], ''))
    elif args.table and not args.raw:
        if len(args.analyses) != 1:
            stderr("Can only request table data on one Analysis at a time.")
            sys.exit(1)
        api_helper(args, route="analyses", supplement="/table")


def references(args):
    api_helper(args, route="references")
