"""
Functions to implement the v0 One Codex API calls.
"""
import json
import os
import requests
from requests_toolbelt import MultipartEncoder, MultipartEncoderMonitor
import sys
from multiprocessing import Lock, Value
from threading import BoundedSemaphore, Thread
import urlparse
from onecodex import version


# Config
if os.environ.get("ONE_CODEX_API_BASE") is not None:
    BASE_API = os.environ.get("ONE_CODEX_API_BASE")
    print "ALL REQUESTS GOING THROUGH: %s" % BASE_API
else:
    BASE_API = "https://beta.onecodex.com/api/v0/"

BASE_URL = urlparse.urlparse(BASE_API)
BASE_URL = BASE_URL._replace(path='/').geturl()
DEFAULT_THREADS = 4
CHUNK_SIZE = 8192

BAD_AUTH_MSG = ("\nYour login credentials appear be bad. Try logging out:"
                "\n    onecodex logout"
                "\n"
                "\nAnd then logging back in:"
                "\n    onecodex login"
                "\n")


# Helpers
def pprint(j, args):
    if args.pprint:
        print json.dumps(j, sort_keys=True,
                         indent=4, separators=(',', ': '))
    else:
        print j


def download_file_helper(url, input_path, auth=None):
    r = requests.get(url, stream=True, auth=auth)
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
    print "Successfully downloaded %s to %s" % (original_filename, local_full_path)


# Version checking function
def get_update_message():
    r = requests.post(BASE_API + "check_for_cli_update",
                      data={"version": version.VERSION,
                            "api_version": version.API_VERSION})
    if r.status_code == 200:
        j = r.json()
        if j.get("message"):
            print j["message"]


def upload_callback(monitor, upload_progress_bytes, lock, total_bytes, n_files):
    if upload_progress_bytes.value == -1:
        return
    with lock:
        upload_progress_bytes.value += CHUNK_SIZE  # Chunk size
    if upload_progress_bytes.value == 0:
        progress = 0.0
    else:
        progress = upload_progress_bytes.value / float(total_bytes)
    # print upload_progress_bytes.value, monitor.bytes_read
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
    sys.stdout.write(text)
    sys.stdout.flush()


# Upload functions
def upload(args):
    """
    Note that this doesn't actually use the default API route -- it instead
    posts directly to S3.
    """
    creds = (args.credentials['api_key'], '')

    if args.threads:
        semaphore = BoundedSemaphore(args.max_threads)
        if args.max_threads != DEFAULT_THREADS:
            print "Uploading with up to %d threads." % args.max_threads

    # Get the initially needed routes
    r0 = requests.get(BASE_API + 'presign_upload', auth=creds)
    if r0.status_code == 401:
        print BAD_AUTH_MSG
        sys.exit(1)
    elif r0.status_code != 200:
        print "Failed to get upload signing credentials"
        sys.exit(1)

    j0 = r0.json()
    s3_url = j0['url']
    signing_url = BASE_URL.rstrip("/") + j0['signing_url']
    callback_url = BASE_URL.rstrip("/") + j0['callback_url']

    upload_threads = []
    upload_progress_bytes = Value('i', 0)
    upload_progress_lock = Lock()
    total_bytes = sum([os.path.getsize(f) for f in args.file])
    total_files = Value('i', len(args.file))
    for f in args.file:
        if args.threads and len(args.file) > 1:  # parallel uploads
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
        print "Failed to get upload signing credentials"
        sys.exit(1)

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
    r2 = requests.post(s3_url, data=m, headers={"Content-Type": m.content_type})
    if r2.status_code != 201:
        print "Upload failed. Please contact help@onecodex.com for assistance."
        sys.exit(1)

    # Finally, issue a callback
    r3 = requests.post(callback_url, auth=creds, data={
        "location": r2.headers['location'],
        "size": os.path.getsize(f)
    })
    if r3.status_code == 200:
        if upload_progress_bytes.value > 0:  # == -1 upon completion
            print "\r"
        print "Successfully uploaded: %s" % f
        with upload_progress_lock:
            total_files.value -= 1
    else:
        print "Failed to upload: %s" % f
        sys.exit(1)

    if semaphore is not None:
        semaphore.release()


# Helper for /route/UUID pattern
def api_helper(args, route, supplement=""):
    creds = (args.credentials['api_key'], '')
    if not getattr(args, route):
        r = requests.get(BASE_API + route + supplement,
                         auth=creds)
        j = r.json()
        pprint(j, args)
    else:
        for uuid in getattr(args, route):
            r = requests.get(BASE_API + route + "/" + uuid + supplement,
                             auth=creds)
            j = r.json()
            pprint(j, args)

    if r.status_code == 401:
        print BAD_AUTH_MSG
        sys.exit(1)


def samples(args):
    api_helper(args, route="samples")


def analyses(args):
    if not args.raw and not args.table:
        api_helper(args, route="analyses")
    elif args.raw and args.table:
        print "Can only request raw or table data at the same time."
        sys.exit(1)
    elif args.raw and not args.table:
        if len(args.analyses) != 1:
            print "Can only request raw data on one Analysis at a time."
            sys.exit(1)
        download_file_helper(BASE_API + "analyses/" + args.analyses[0] + "/raw",
                             input_path=args.raw,
                             auth=(args.credentials['api_key'], ''))
    elif args.table and not args.raw:
        if len(args.analyses) != 1:
            print "Can only request table data on one Analysis at a time."
            sys.exit(1)
        api_helper(args, route="analyses", supplement="/table")


def references(args):
    api_helper(args, route="references")
