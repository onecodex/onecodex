"""
Functions for connecting to the One Codex server; these should be rolled out
into the onecodex python library at some point for use across CLI and GUI clients
"""
from __future__ import print_function, division

import bz2
import click
from collections import deque, OrderedDict
import gzip
from multiprocessing import Value
import os
import re
import requests
from requests_toolbelt import MultipartEncoder
from threading import BoundedSemaphore, Thread
from unidecode import unidecode
import warnings

from onecodex.exceptions import OneCodexException, UploadException


DEFAULT_UPLOAD_THREADS = 4


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
        ret = b''.join(ret_list)
        self._size -= len(ret)
        return ret

    def flush(self):
        pass

    def close(self):
        self.closed = True


class FASTXInterleave(object):
    """Wrapper around two `file` objects that decompresses gzip or bz2, where applicable, and
    interleaves the two files either two or four lines at a time. Yields uncompressed data.

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

    def __init__(self, file_path, file_size, file_format='fastq', progressbar=None):
        if file_path[0].endswith('.gz') or file_path[1].endswith('.gz'):
            self._fp_left = gzip.GzipFile(file_path[0], mode='rb')
            self._fp_right = gzip.GzipFile(file_path[1], mode='rb')
        elif file_path[0].endswith('.bz2') or file_path[1].endswith('.bz2'):
            self._fp_left = bz2.BZ2File(file_path[0], mode='rb')
            self._fp_right = bz2.BZ2File(file_path[1], mode='rb')
        else:
            self._fp_left = open(file_path[0], mode='rb')
            self._fp_right = open(file_path[1], mode='rb')

        if file_format == 'fasta':
            self._lines_per_record = 2
        elif file_format == 'fastq':
            self._lines_per_record = 4
        else:
            raise OneCodexException('file_format must be one of: fastq, fasta')

        self._tell = 0
        self._max_tell = 0
        self._fsize = file_size
        self._buf = Buffer()

        self.progressbar = progressbar
        self.mime_type = 'text/plain'

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
            # only update bar when we move further in the file. we might have done a seek(0) to
            # retry the upload, in which case we don't want to double-count our progress
            if self._tell > self._max_tell:
                self.progressbar.update(self._tell - self._max_tell)

        if self._tell > self._max_tell:
            self._max_tell = self._tell

        return bytes_read

    def seek(self, loc):
        """Called if upload fails and must be retried."""
        assert loc == 0

        self._fp_left.seek(loc)
        self._fp_right.seek(loc)
        self._tell = loc
        self._buf = Buffer()

    def close(self):
        self._fp_left.close()
        self._fp_right.close()


class FASTXPassthru(object):
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
        self._fp = open(file_path, mode='rb')
        self._fsize = file_size
        self._max_tell = 0
        self.progressbar = progressbar

        _, ext = os.path.splitext(file_path)

        if ext in {'.gz', '.gzip'}:
            self.mime_type = 'application/x-gzip'
        elif ext in {'.bz', '.bz2', '.bzip', '.bzip2'}:
            self.mime_type = 'application/x-bzip2'
        else:
            self.mime_type = 'text/plain'

    def read(self, size=-1):
        bytes_read = self._fp.read(size)

        if self.progressbar:
            # only update bar when we move further in the file. we might have done a seek(0) to
            # retry the upload, in which case we don't want to double-count our progress
            if self._fp.tell() > self._max_tell:
                self.progressbar.update(self._fp.tell() - self._max_tell)

        if self._fp.tell() > self._max_tell:
            self._max_tell = self._fp.tell()

        return bytes_read

    @property
    def len(self):
        """Size of data left to be read."""
        return self._fsize - self._fp.tell()

    def seek(self, loc):
        """Called if upload fails and must be retried."""
        assert loc == 0

        self._fp.seek(loc)

    def close(self):
        self._fp.close()


def interleaved_filename(file_path):
    """Return filename used to represent a set of paired-end files. Assumes Illumina-style naming
    conventions where each file has _R1_ or _R2_ in its name."""
    if not isinstance(file_path, tuple):
        raise OneCodexException('Cannot get the interleaved filename without a tuple.')
    if re.match('.*[._][Rr][12][_.].*', file_path[0]):
        return re.sub('[._][Rr][12]', '', file_path[0])
    else:
        warnings.warn('Paired-end filenames do not match--are you sure they are correct?')
        return file_path[0]


def _file_size(file_path, uncompressed=False):
    """Return size of a single file, compressed or uncompressed"""
    _, ext = os.path.splitext(file_path)

    if uncompressed:
        if ext in {'.gz', '.gzip'}:
            with gzip.GzipFile(file_path, mode='rb') as fp:
                try:
                    fp.seek(0, os.SEEK_END)
                    return fp.tell()
                except ValueError:
                    # on python2, cannot seek from end and must instead read to end
                    fp.seek(0)
                    while len(fp.read(8192)) != 0:
                        pass
                    return fp.tell()
        elif ext in {'.bz', '.bz2', '.bzip', '.bzip2'}:
            with bz2.BZ2File(file_path, mode='rb') as fp:
                fp.seek(0, os.SEEK_END)
                return fp.tell()
        else:
            return os.path.getsize(file_path)
    else:
        return os.path.getsize(file_path)


def _file_stats(file_path):
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
    else:
        file_size = _file_size(file_path, uncompressed=False)

    new_filename, ext = os.path.splitext(os.path.basename(file_path))

    if ext in {'.gz', '.gzip', '.bz', '.bz2', '.bzip'}:
        new_filename, ext = os.path.splitext(new_filename)

    final_filename = new_filename + ext

    if ext in {'.fa', '.fna', '.fasta'}:
        file_format = 'fasta'
    elif ext in {'.fq', '.fastq'}:
        file_format = 'fastq'
    else:
        raise UploadException(
            '{}: extension must be one of .fa, .fna, .fasta, .fq, .fastq'.format(final_filename)
        )

    if file_size == 0:
        raise UploadException('{}: empty files can not be uploaded'.format(final_filename))

    return final_filename, file_size, file_format


# this lets us turn off the click progressbar context manager and is python2 compatible
# https://stackoverflow.com/questions/45187286/how-do-i-write-a-null-no-op-contextmanager-in-python
class FakeProgressBar(object):
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


def process_api_error(resp, state=None):
    """Raise an exception with a pretty message in various states of upload"""
    error_code = resp.status_code
    error_json = resp.json()

    if error_code == 402:
        error_message = ('Please add a payment method to upload more samples. If you continue to '
                         'experience problems, contact us at help@onecodex.com for assistance.')
    else:
        if 'msg' in error_json:
            error_message = error_json['msg'].rstrip('.')
        elif 'message' in error_json:
            error_message = error_json['message'].rstrip('.')
        else:
            error_message = None

        if state == 'init' and not error_message:
            error_message = ('Could not initialize upload. Are you logged in? If this problem '
                             'continues, please contact help@onecodex.com for assistance.')
        elif state == 'upload' and not error_message:
            error_message = ('File could not be uploaded. If this problem continues, please contact '
                             'help@onecodex.com for assistance.')
        elif state == 'callback' and not error_message:
            error_message = ('Callback could not be completed. If this problem continues, please '
                             'contact help@onecodex.com for assistance.')

    if error_message is None:
        error_message = 'Upload failed. Please contact help@onecodex.com for assistance.'

    raise UploadException(error_message)


def upload(files, session, samples_resource, threads=DEFAULT_UPLOAD_THREADS, log=None,
           metadata=None, tags=None, project=None, coerce_ascii=None, progressbar=False):
    """Uploads multiple files to the One Codex server via either fastx-proxy or directly to S3.

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
    threads : `integer`, optional
        Number of concurrent uploads. May provide a speedup.
    log : `logging.Logger`, optional
        Used to write status messages to a file or terminal.
    metadata : `dict`, optional
    tags : `list`, optional
    project : `string`, optional
        UUID of project to associate this sample with.
    coerce_ascii : `bool`, optional
        If true, rename unicode filenames to ASCII and issue warning.
    progressbar : `bool`, optional
        If true, display a progress bar using Click.

    Returns
    -------
    `list` of `string`s containing sample UUIDs of newly uploaded files.
    """
    filenames = []
    file_sizes = []
    file_formats = []

    for file_path in files:
        new_filename, file_size, file_format = _file_stats(file_path)
        filenames.append(new_filename)
        file_sizes.append(file_size)
        file_formats.append(file_format)

    # if filename cannot be represented as ascii, raise and suggest renaming
    for idx, fname in enumerate(filenames):
        ascii_fname = unidecode(fname)

        if fname != ascii_fname:
            if coerce_ascii:
                if log:
                    log.warning('Renaming {} to {}, must be ASCII\n'.format(
                        fname.encode('utf-8'), ascii_fname
                    ))
                filenames[idx] = ascii_fname
            else:
                raise OneCodexException('Filenames must be ascii. Try using --coerce-ascii')

    # threaded upload setup
    upload_threads = []
    uploaded_uuids = []

    if threads > 1:
        import ctypes
        thread_error = Value(ctypes.c_wchar_p, '')
        semaphore = BoundedSemaphore(threads)

        def threaded_upload(*args):
            def _wrapped(*wrapped_args):
                semaphore.acquire()

                try:
                    file_uuid = upload_fileobj(*wrapped_args[:-1])
                    uploaded_uuids.append(file_uuid)
                except Exception as e:
                    # handle inside the thread to prevent the exception message from leaking out
                    wrapped_args[-1].value = '{}'.format(e)

                    if log:
                        log.error('{}: {}'.format(wrapped_args[1], e))

                semaphore.release()

            # the thread error message must be the last parameter
            thread = Thread(target=_wrapped, args=args + (thread_error, ))
            thread.daemon = True
            thread.start()
            upload_threads.append(thread)

    # file_path is the path to the file on this disk. filename is what we'll call the file in the
    # mainline database. file_size is the sum of both files in a pair, or the size of an unpaired
    # file. if paired, file_size is the uncompressed size. if unpaired, file_size is the actual
    # size on disk. unpaired files are uploaded as-is. paired files are decompressed, interleaved,
    # and uploaded as uncompressed data.
    if progressbar:
        bar_context = click.progressbar(length=sum(file_sizes), label='Uploading... ')
    else:
        bar_context = FakeProgressBar()

    with bar_context as bar:
        for file_path, filename, file_size, file_format in zip(files, filenames, file_sizes, file_formats):
            if isinstance(file_path, tuple):
                file_obj = FASTXInterleave(file_path, file_size, file_format, bar)
            else:
                file_obj = FASTXPassthru(file_path, file_size, bar)

            if threads > 1:
                threaded_upload(
                    file_obj, filename, file_size, session, samples_resource, log, metadata, tags, project
                )
            else:
                try:
                    file_uuid = upload_fileobj(
                        file_obj, filename, file_size, session, samples_resource, log, metadata, tags, project
                    )
                    uploaded_uuids.append(file_uuid)
                except UploadException as e:
                    if log:
                        log.error('{}: {}'.format(filename, e))

        if threads > 1:
            # we need to do this funky wait loop to ensure threads get killed by ctrl-c
            while True:
                for thread in upload_threads:
                    # hopefully no one has a file that takes longer than a week to upload
                    thread.join(604800)

                if all(not thread.is_alive() for thread in upload_threads):
                    break

        bar.finish()

    return uploaded_uuids


def upload_fileobj(file_obj, filename, file_size, session, samples_resource, log=None, metadata=None,
                   tags=None, project=None):
    """Uploads a single file-like object to the One Codex server either via fastx-proxy or directly
    to S3.

    Parameters
    ----------
    file_obj : `FASTXInterleave`, `FASTXPassthru`, or a file-like object
        A wrapper around a pair of fastx files (`FASTXInterleave`) or a single fastx file. In the
        case of paired files, they will be interleaved and uploaded uncompressed. In the case of a
        single file, it will simply be passed through (`FASTXPassthru`) to One Codex, compressed
        or otherwise. If a file-like object is given, its mime-type will be sent as 'text/plain'.
    filename : `string`
        The filename you wish to associate this fastx file with at One Codex.
    file_size : `integer`
        Accurate size of file to be uploaded, in bytes.
    session : `requests.Session`
        Connection to One Codex API.
    samples_resource : `onecodex.models.Samples`
        Wrapped potion-client object exposing `init_upload` and `confirm_upload` routes to mainline.
    log : `logging.Logger`, optional
        Used to write status messages to a file or terminal.
    metadata : `dict`, optional
    tags : `list`, optional
    project : `string`, optional
        UUID of project to associate this sample with.

    Returns
    -------
    `string` containing sample UUID of newly uploaded file.
    """
    upload_args = {
        'filename': filename,
        'size': file_size,
        'upload_type': 'standard'  # this is multipart form data
    }

    if metadata:
        upload_args['metadata'] = metadata

    if tags:
        upload_args['tags'] = tags

    if project:
        upload_args['project'] = project.id

    try:
        upload_info = samples_resource.init_upload(upload_args)
    except requests.exceptions.HTTPError as e:
        process_api_error(e.response, state='init')
    except requests.exceptions.ConnectionError:
        raise UploadException(
            'The command line client is experiencing connectivity issues and '
            'cannot complete the upload of %s at this time. Please try again '
            'later. If the problem persists, contact us at help@onecodex.com '
            'or assistance.' % filename
        )

    upload_url = upload_info['upload_url']

    # need an OrderedDict to preserve order for S3--does this matter?
    multipart_fields = OrderedDict()

    for k, v in upload_info['additional_fields'].items():
        multipart_fields[str(k)] = str(v)

    # this attribute is only in FASTXInterleave and FASTXPassthru
    mime_type = getattr(file_obj, 'mime_type', 'text/plain')

    multipart_fields['file'] = (filename, file_obj, mime_type)
    encoder = MultipartEncoder(multipart_fields)
    content_type = encoder.content_type

    # try to upload the file, retrying as necessary
    max_retries = 3
    n_retries = 0

    while n_retries < max_retries:
        try:
            upload_request = session.post(
                upload_url, data=encoder, headers={'Content-Type': content_type}, auth={}
            )

            if upload_request.status_code not in [200, 201]:
                # if using proxy, special route can provide more info on e.g., validation issues
                if multipart_fields.get('sample_id'):
                    error_url = '/'.join(upload_url.split('/')[:-1]) + '/errors'

                    try:
                        e_resp = session.post(error_url, json={'sample_id': multipart_fields.get('sample_id')})

                        if e_resp.status_code == 200:
                            process_api_error(e_resp, state='upload')
                        else:
                            process_api_error(upload_request, state='upload')
                    except requests.exceptions.RequestException:
                        process_api_error(upload_request, state='upload')

            file_obj.close()
            break
        except requests.exceptions.ConnectionError:
            n_retries += 1
            file_obj.seek(0)  # reset file_obj back to start

            if n_retries == max_retries:
                raise UploadException(
                    'The command line client is experiencing connectivity issues and '
                    'cannot complete the upload of %s at this time. Please try again '
                    'later. If the problem persists, contact us at help@onecodex.com '
                    'or assistance.' % filename
                )
            else:
                if log:
                    log.error(
                        '{}: Connectivity issue. Retry {} of {}.'.format(filename, n_retries, max_retries)
                    )

    # finally, issue a callback
    try:
        if not multipart_fields.get('callback_url'):
            samples_resource.confirm_upload({
                'sample_id': upload_info['sample_id'],
                'upload_type': 'standard'
            })
    except requests.exceptions.HTTPError as e:
        process_api_error(e.response, state='callback')

    if log:
        log.info('{}: finished as sample {}'.format(filename, upload_info['sample_id']))

    return upload_info['sample_id']
