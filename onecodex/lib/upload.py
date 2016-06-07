"""
Functions for connecting to the One Codex server; these should be rolled out
into the onecodex python library at some point for use across CLI and GUI clients
"""
from __future__ import print_function

from math import floor
import os

import boto3
from boto3.s3.transfer import S3Transfer
from boto3.exceptions import S3UploadFailedError


class UploadException(Exception):
    """
    A exception for when things go wrong with uploading
    """
    pass


def upload_file(filename, session, server_url, progress_callback=None, n_callbacks=400):
    """
    Uploads a file to the One Codex server (and handles files >5Gb)

    Takes an optional callback that it calls with a number from 0 to 1 as the
    upload progresses.
    """
    # first check with the one codex server to get upload parameters
    req = session.get(server_url + 'api/v1/samples/init_multipart_upload')
    if req.status_code != 200:
        raise UploadException('Could not initiate upload with One Codex server')

    upload_params = req.json()
    callback_url = server_url.rstrip('/') + upload_params['callback_url']

    access_key = upload_params['upload_aws_access_key_id']
    secret_key = upload_params['upload_aws_secret_access_key']

    # set up a progress tracker to simplify the callback
    if progress_callback is not None:
        class Progress(object):
            """
            Wrapper for progress callbacks
            """
            def __init__(self, callback, file_size):
                self.file_size = float(file_size)
                self.transferred = 0
                self.callback = callback
                self.step_size = n_callbacks  # e.g. call callback in 1000 or less intervals

            def __call__(self, bytes_seen):
                per_prev_done = self.transferred / self.file_size
                self.transferred += bytes_seen
                per_done = self.transferred / self.file_size
                if floor(self.step_size * per_prev_done) != floor(self.step_size * per_done):
                    self.callback(filename, per_done)
        progress_tracker = Progress(progress_callback, os.path.getsize(filename))
    else:
        progress_tracker = None

    # actually do the upload
    client = boto3.client('s3', aws_access_key_id=access_key, aws_secret_access_key=secret_key)
    # TODO: this automatically uses 10 threads, but we'd probably like it to be configurable
    transfer = S3Transfer(client)
    try:
        transfer.upload_file(filename, upload_params['s3_bucket'], upload_params['file_id'],
                             extra_args={'ServerSideEncryption': 'AES256'},
                             callback=progress_tracker)
    except S3UploadFailedError:
        raise UploadException('Upload has failed. Please contact help@onecodex.com '
                              'if you experience further issues')

    # return completed status to the one codex server
    s3_path = 's3://{}/{}'.format(upload_params['s3_bucket'], upload_params['file_id'])
    req = session.post(callback_url,
                       json={'s3_path': s3_path, 'filename': os.path.basename(filename)})

    if req.status_code != 200:
        raise UploadException('Upload confirmation has failed. Please contact help@onecodex.com '
                              'if you experience further issues')
