import bz2
from collections import deque
import gzip
from io import BytesIO
import os
import re
import string
import warnings

from onecodex.exceptions import ValidationError, ValidationWarning

GZIP_COMPRESSION_LEVEL = 5


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


class GzipBuffer(object):
    def __init__(self):
        self._buf = Buffer()
        self._gzip = gzip.GzipFile(None, mode='wb', fileobj=self._buf,
                                   compresslevel=GZIP_COMPRESSION_LEVEL)
        self._reads_buffer = Buffer()
        self.MAX_READS_BUFFER_SIZE = 1024 * 256  # 256kb
        self.closed = False

    def __len__(self):
        return len(self._buf)

    def write(self, s):
        self._reads_buffer.write(s)
        if len(self._reads_buffer) >= self.MAX_READS_BUFFER_SIZE:
            self.flush()

    def read(self, size=-1):
        return self._buf.read(size)

    def flush(self):
        self._gzip.write(self._reads_buffer.read())

    def close(self):
        if len(self._reads_buffer) > 0:
            self.flush()
        self._gzip.close()
        self.closed = True


# this checks and translates all valid IUPAC nucleotide codes into the core 4+n (ACGTN)
OTHER_BASES = re.compile(b'[BDHIKMRSUVWXYbdhikmrsuvwxy]')
if hasattr(bytes, 'maketrans'):
    OTHER_BASE_TRANS = bytes.maketrans(b'BDHIKMRSUVWXYbdhikmrsuvwxy',
                                       b'NNNNNNNNTNNNNnnnnnnnntnnnn')
else:
    OTHER_BASE_TRANS = string.maketrans(b'BDHIKMRSUVWXYbdhikmrsuvwxy',
                                        b'NNNNNNNNTNNNNnnnnnnnntnnnn')


class FASTXNuclIterator(object):
    def __init__(self, file_obj, allow_iupac=False, check_filename=True, as_raw=False,
                 validate=True):
        if hasattr(file_obj, 'name'):
            self.name = file_obj.name
        else:
            self.name = 'File'

        self._set_file_obj(file_obj, check_filename=check_filename)
        self.unchecked_buffer = b''
        self.seq_reader = self._generate_seq_reader(False)
        self.allow_iupac = allow_iupac
        self.validate = validate

        if self.allow_iupac:
            self.valid_bases = re.compile(b'[^ABCDGHIKMNRSTUVWXYabcdghikmnrstuvwxy\s]')
            self.valid_bases_match = re.compile(b'^[ABCDGHIKMNRSTUVWXYabcdghikmnrstuvwxy\s]*$')
        else:
            self.valid_bases = re.compile(b'[^ACGTNacgtn\s]')
            self.valid_bases_match = re.compile(b'^[ACGTNacgtn\s]*$')
        self.as_raw = as_raw

        self._set_total_size()
        self.processed_size = self.file_obj.tell()
        self.warnings = set()

    def _set_file_obj(self, file_obj, check_filename=True):
        """
        Transparently decompress files and determine what kind of file they are (FASTA/Q).
        """
        if not hasattr(file_obj, 'name'):
            # can't do the checks if there's not filename
            check_filename = False

        # detect if gzipped/bzipped and uncompress transparently
        start = file_obj.read(1)
        if start == b'\x1f':
            if check_filename and not file_obj.name.endswith(('.gz', '.gzip')):
                raise ValidationError('{} is gzipped, but lacks a ".gz" ending'.format(self.name))
            file_obj.seek(0)
            file_obj = gzip.GzipFile(fileobj=file_obj)
            start = file_obj.read(1)
        elif start == b'\x42' and hasattr(bz2, 'open'):
            if check_filename and not file_obj.name.endswith(('.bz2', '.bz', '.bzip')):
                raise ValidationError('{} is bzipped, but lacks a ".bz2" ending'.format(self.name))
            # we can only read BZ2 files in python 3.3 and above
            file_obj.seek(0)
            file_obj = bz2.open(file_obj)
            start = file_obj.read(1)

        # determine if a FASTQ or a FASTA
        if start == b'>':
            self.file_type = 'FASTA'
            if check_filename and not ('.fa' in file_obj.name or
                                       '.fna' in file_obj.name or
                                       '.fasta' in file_obj.name):
                raise ValidationError('{} is FASTA, but lacks a ".fa" ending'.format(self.name))
        elif start == b'@':
            self.file_type = 'FASTQ'
            if check_filename and not ('.fq' in file_obj.name or
                                       '.fastq' in file_obj.name):
                raise ValidationError('{} is FASTQ, but lacks a ".fq" ending'.format(self.name))
        else:
            raise ValidationError('{} is not valid FASTX'.format(self.name))

        self.file_obj = file_obj

    def _set_total_size(self):
        if isinstance(self.file_obj, BytesIO):
            self.file_obj.seek(0)
            self.total_size = len(self.file_obj.read())
            self.file_obj.seek(1)
        else:
            try:
                self.total_size = os.fstat(self.file_obj.fileno()).st_size
                if self.total_size < 70:
                    raise ValidationError('{} is too small to be analyzed: {} bytes'.format(
                        self.name, self.total_size
                    ))
            except IOError:
                pass

        # Set the buffer size, 16MB by default for files >32MB
        if self.total_size >= (1024 * 1024 * 32):
            self.buffer_read_size = 1024 * 1024 * 16  # 16MB
        else:
            self.buffer_read_size = 1024 * 16  # 16KB small chunk

    def _generate_seq_reader(self, last=False):
        # the last record doesn't have a @/> on the next line so we omit that
        # if the "last" flag is passed (to allow reading the last record)
        if self.file_type == 'FASTA':
            seq_reader = re.compile(b"""
                (?P<id>[^\\n]+)\\n  # the identifier line
                (?P<seq>[^>]+)  # the sequence
            """ + (b'' if last else b'(?:\\n>)'), re.VERBOSE)
        elif self.file_type == 'FASTQ':
            seq_reader = re.compile(b"""
                (?P<id>[^\\n]+)\\n
                (?P<seq>[^\\n]+)\\n
                \+(?P<id2>[^\\n]*)\\n
                (?P<qual>[^\\n]+)
            """ + (b'' if last else b'(?:\\n@)'), re.DOTALL + re.VERBOSE)
        return seq_reader

    def _warn_once(self, message):
        if message in self.warnings:
            return
        warnings.warn(message, ValidationWarning)
        self.warnings.add(message)

    def _validate_record(self, rec):
        # TODO: if there are quality scores, make sure they're in range
        # FIXME: fail if reads aren't interleaved and an override flag isn't passed?
        seq_id, seq, seq_id2, qual = rec['id'], rec['seq'], rec.get('id2', b''), rec.get('qual')
        if not self.validate:
            return seq_id, seq, seq_id2, qual

        if b'\t' in seq_id or b'\t' in seq_id2:
            self._warn_once('{} can not have tabs in headers; autoreplacing'.format(self.name))
            seq_id = seq_id.replace('\t', '|')

        # Match then search is ~5-10% faster than just searching
        if not bool(self.valid_bases_match.match(seq)):
            chars = ','.join(set(self.valid_bases.findall(seq)))
            raise ValidationError('{} contains non-nucleic acid characters: {}'.format(self.name,
                                                                                       chars))

        # Only search for OTHER_BASES if we're allowing them above in the first place
        if self.allow_iupac and OTHER_BASES.search(seq) is not None:
            self._warn_once('Translating other bases in {} (X->N,U->T)'.format(self.name))
            seq = seq.translate(OTHER_BASE_TRANS)

        return seq_id, seq, seq_id2, qual

    def __iter__(self):
        eof = False
        while not eof:
            new_data = self.file_obj.read(self.buffer_read_size)
            # if we're at the end of the file
            if len(new_data) == 0:
                # switch to a different regex to parse without a next record
                eof = True
                self.seq_reader = self._generate_seq_reader(True)
                # automatically remove newlines from the end of the file (they get added back in
                # by the formatting operation below, but otherwise they mess up the regex and you
                # end up with two terminating \n's)
                self.unchecked_buffer = self.unchecked_buffer.rstrip(b'\n')
            else:
                self.unchecked_buffer += new_data

            end = 0
            while True:
                match = self.seq_reader.match(self.unchecked_buffer, end)
                if match is None:
                    break
                rec = match.groupdict()
                seq_id, seq, seq_id2, qual = self._validate_record(rec)
                if self.as_raw:
                    yield (seq_id, seq, qual)
                elif self.file_type == 'FASTA':
                    yield b'>' + seq_id + b'\n' + seq + b'\n'
                elif self.file_type == 'FASTQ':
                    yield (b'@' + seq_id + b'\n' + seq +
                           b'\n+' + seq_id2 + b'\n' + qual + b'\n')
                end = match.end()

            if hasattr(self.file_obj, 'fileobj'):
                # for gzip files, get the amount read of the gzipped file (which is wrapped inside)
                self.processed_size = self.file_obj.fileobj.tell()
            else:
                self.processed_size = self.file_obj.tell()
            self.unchecked_buffer = self.unchecked_buffer[end:]

    @property
    def bytes_left(self):
        if self.total_size is not None:
            return self.total_size - self.processed_size

    def close(self):
        # did we read everything?
        if self.bytes_left != 0:
            raise ValidationError('Failed to properly read file: {}/{} bytes unread.'.format(
                self.bytes_left, self.total_size))

        # actually close the file
        self.file_obj.close()


class FASTXTranslator(object):
    def __init__(self, file_obj, pair=None, recompress=True, progress_callback=None,
                 total=None, **kwargs):
        # detect if gzipped/bzipped and uncompress transparently
        self.reads = FASTXNuclIterator(file_obj, **kwargs)
        self.reads_iter = iter(self.reads)
        if pair is not None:
            self.reads_pair = FASTXNuclIterator(pair)
            self.reads_pair_iter = iter(self.reads_pair)
            if self.reads.file_type != self.reads_pair.file_type:
                raise ValidationError('Paired read files are different types (FASTA/FASTQ)')
        else:
            self.reads_pair = None
            self.reads_pair_iter = None

        if recompress:
            self.checked_buffer = GzipBuffer()
        else:
            self.checked_buffer = Buffer()

        self.progress_callback = progress_callback
        self.total = total
        self.total_written = 0

        # save in case we need to reset later
        self._saved_args = kwargs.copy()
        self._saved_args.update({'recompress': recompress, 'progress_callback': progress_callback})

    def read(self, n=-1):
        if self.reads_pair is None:
            while len(self.checked_buffer) < n or n < 0:
                try:
                    record = next(self.reads_iter)
                except StopIteration:
                    record = None

                if record is not None:
                    self.checked_buffer.write(record)
                elif record is None:
                    self.checked_buffer.close()
                    break

                if self.progress_callback is not None:
                    self.progress_callback(self.reads.name, self.reads.processed_size,
                                           validation=(not self.reads.validate))
        else:
            while len(self.checked_buffer) < n or n < 0:
                try:
                    record = next(self.reads_iter)
                except StopIteration:
                    record = None
                try:
                    record_pair = next(self.reads_pair_iter)
                except StopIteration:
                    record_pair = None

                if record is not None and record_pair is not None:
                    self.checked_buffer.write(record)
                    self.checked_buffer.write(record_pair)
                elif record is None and record_pair is None:
                    self.checked_buffer.close()
                    break
                else:
                    raise ValidationError("Paired read files do not have the "
                                          "same number of records")

                if self.progress_callback is not None:
                    bytes_uploaded = self.reads.processed_size + self.reads_pair.processed_size
                    self.progress_callback(self.reads.name, bytes_uploaded,
                                           validation=(not self.reads.validate))

        bytes_reads = self.checked_buffer.read(n)
        self.total_written += len(bytes_reads)
        return bytes_reads

    def readall(self):
        return self.read()

    @property
    def len(self):
        # Properly ensure requests_toolbelt reads the entirety of the
        # file *plus* the remaining buffer object
        if self.total is None:
            self.reads.validate = False
            if self.reads_pair:
                self.reads_pair.validate = False
            while len(self.read(8192)) != 0:
                pass
            self.total = self.total_written
            self.seek(0)
        return self.total - self.total_written

    def __len__(self):
        return self.len

    def seek(self, loc):
        assert loc == 0  # we can only rewind all the way
        reads = self.reads.file_obj
        reads.seek(0)
        if self.reads_pair:
            pair = self.reads_pair.file_obj
            pair.seek(0)
        else:
            pair = None
        self.__init__(reads, pair, total=self.total, **self._saved_args)

    def write(self, b):
        raise NotImplementedError()

    def close(self):
        assert len(self.checked_buffer) == 0
        self.reads.close()
        if self.reads_pair is not None:
            self.reads_pair.close()
