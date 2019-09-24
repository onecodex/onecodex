import os
import re
import logging
from unidecode import unidecode

from onecodex.exceptions import OneCodexException, UploadException

R1_FILENAME_RE = re.compile(".*[._][Rr]?[1][_.].*")
R2_FILENAME_RE = re.compile(".*[._][Rr]?[2][_.].*")
log = logging.getLogger("onecodex")


def _check_for_ascii_filename(filename, coerce_ascii):
    """Check that the filename is ASCII.

    If it isn't, convert it to ASCII & return it if the ascii flag
    has been set otherwise raise an exception.
    """
    try:
        # python2
        ascii_fname = unidecode(unicode(filename))
    except NameError:
        ascii_fname = unidecode(filename)

    if filename != ascii_fname:
        if coerce_ascii:
            # TODO: Consider warnings.warn here instead
            log.warning(
                "Renaming {} to {}, must be ASCII\n".format(filename.encode("utf-8"), ascii_fname)
            )
            filename = ascii_fname
        else:
            raise OneCodexException("Filenames must be ascii. Try using --coerce-ascii")
    return filename


def get_fastx_format(file_path):
    """Return format of given file: fasta or fastq.

    Assumes Illumina-style naming conventions where each file has _R1_ or _R2_ in its name.
    If the file is not fasta or fastq, raises an exception
    """
    new_filename, ext = os.path.splitext(os.path.basename(file_path))

    if ext in {".gz", ".gzip", ".bz", ".bz2", ".bzip"}:
        new_filename, ext = os.path.splitext(new_filename)

    if ext in {".fa", ".fna", ".fasta"}:
        return "fasta"
    elif ext in {".fq", ".fastq"}:
        return "fastq"
    else:
        raise UploadException(
            "{}: extension must be one of .fa, .fna, .fasta, .fq, .fastq".format(file_path)
        )


class FilePassthru(object):
    """Wrapper around `file` object that updates a progress bar and guesses mime-type.

    Parameters
    ----------
    file_path : `string`
        Path to file.
    progressbar : `click.progressbar`, optional
        The progress bar to update.
    """

    def __init__(self, file_path, progressbar=None):
        self._fp = open(file_path, mode="rb")
        self._fsize = os.path.getsize(file_path)

        self.progressbar = progressbar

        _, ext = os.path.splitext(file_path)
        self.filename = os.path.basename(file_path)

        if self._fsize == 0:
            raise UploadException("{}: empty files can not be uploaded".format(self.filename))

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

    def size(self):
        return self._fsize

    @property
    def len(self):
        """Size of data left to be read."""
        return self._fsize - self._fp.tell()

    def seek(self, loc):
        """Seek to a position in the file.

        Notes
        -----
        This is called if an upload fails and must be retried.
        """
        assert loc == 0

        # rewind progress bar
        if self.progressbar:
            self.progressbar.update(-self._fp.tell())

        self._fp.seek(loc)

    def close(self):
        self._fp.close()

    def enforce_ascii_filename(self, coerce_ascii):
        """Update the filename to be ASCII. Raises an exception if `coerce_ascii` is `False` and the filename is not ASCII."""
        self.filename = _check_for_ascii_filename(self.filename, coerce_ascii)


class PairedEndFiles(object):
    def __init__(self, files, progressbar=None):
        if len(files) != 2:
            raise OneCodexException("Paired files uploading can only take 2 files")

        for f in files:
            if get_fastx_format(f) != "fastq":
                raise OneCodexException("Interleaving FASTA files is currently unsupported")

        if R1_FILENAME_RE.match(files[0]) and R2_FILENAME_RE.match(files[1]):
            file1 = files[0]
            file2 = files[1]
        elif R2_FILENAME_RE.match(files[0]) and R1_FILENAME_RE.match(files[1]):
            file1 = files[1]
            file2 = files[0]
        else:
            raise OneCodexException("Paired files need to have _R1/_1 and _R2/_2 in their name")

        self.r1 = FilePassthru(file1, progressbar)
        self.r2 = FilePassthru(file2, progressbar)

    def enforce_ascii_filename(self, coerce_ascii):
        self.r1.enforce_ascii_filename(coerce_ascii)
        self.r2.enforce_ascii_filename(coerce_ascii)


def get_file_wrapper(file, coerce_ascii, bar):
    """Take a str or tuple (str) and return the corresponding file wrapper object.

    If there is more than one file, it must be a paired end uploads and the filenames will be validated.
    """
    if isinstance(file, tuple):
        fobj = PairedEndFiles(file, bar)
        fobj.enforce_ascii_filename(coerce_ascii)
        return fobj

    fobj = FilePassthru(file, bar)
    fobj.enforce_ascii_filename(coerce_ascii)
    return fobj
