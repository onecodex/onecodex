"""
Functions to implement the v0 One Codex API calls.
"""
import json
import os
import requests
import sys
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


# Helpers
def pprint(j, args):
    if args.pprint:
        print json.dumps(j, sort_keys=True,
                         indent=4, separators=(',', ': '))
    else:
        print j


# Version checking function
def get_update_message():
    r = requests.post(BASE_API + "check_for_cli_update",
                      data={"version": version.VERSION,
                            "api_version": version.API_VERSION})
    if r.status_code == 200:
        j = r.json()
        if j.get("message"):
            print j["message"]


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
    if r0.status_code != 200:
        print "Failed to get upload signing credentials"
        sys.exit(1)

    j0 = r0.json()
    s3_url = j0['url']
    signing_url = BASE_URL.rstrip("/") + j0['signing_url']
    callback_url = BASE_URL.rstrip("/") + j0['callback_url']

    upload_threads = []
    for f in args.file:
        if args.threads:  # parallel uploads
            # Multi-threaded uploads
            t = Thread(target=upload_helper,
                       args=(f, s3_url, signing_url, callback_url, creds, semaphore))
            upload_threads.append(t)
            t.start()
        else:  # serial uploads
            upload_helper(f, s3_url, signing_url, callback_url, creds)

    if args.threads:
        for ut in upload_threads:
            ut.join()


def upload_helper(f, s3_url, signing_url, callback_url, creds,
                  semaphore=None):
    # First get the signing form data
    if semaphore is not None:
        semaphore.acquire()

    r1 = requests.post(signing_url, data={"filename": f},
                       auth=creds)
    if r1.status_code != 200:
        print "Failed to get upload signing credentials"
        sys.exit(1)

    # Then do a multi-part post directly to S3
    fields = {"file": open(f, mode='rb')}
    fields = dict(r1.json().items())
    r2 = requests.post(s3_url, data=fields, files={'file': open(f, mode='rb')})
    if r2.status_code != 201:
        print "Upload failed. Please contact help@onecodex.com for assistance."
        sys.exit(1)

    # Finally, issue a callback
    r3 = requests.post(callback_url, auth=creds, data={
        "location": r2.headers['location'],
        "size": os.path.getsize(f)
    })
    if r3.status_code == 200:
        print "Successfully uploaded: %s" % f
    else:
        print "Failed to upload: %s" % f
        sys.exit(1)

    if semaphore is not None:
        semaphore.release()


# Helper for /route/UUID pattern
def api_helper(args, route):
    creds = (args.credentials['api_key'], '')
    if not getattr(args, route):
        r = requests.get(BASE_API + route,
                         auth=creds)
        j = r.json()
        pprint(j, args)
    else:
        for uuid in getattr(args, route):
            r = requests.get(BASE_API + route + "/" + uuid,
                             auth=creds)
            j = r.json()
            pprint(j, args)


def samples(args):
    api_helper(args, route="samples")


def analyses(args):
    api_helper(args, route="analyses")


def references(args):
    api_helper(args, route="references")
