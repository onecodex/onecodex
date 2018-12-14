import bz2
import gzip
from io import BytesIO
import pytest
import random
import sys
import warnings


from onecodex.exceptions import ValidationError, ValidationWarning
from onecodex.lib.inline_validator import FASTXNuclIterator, FASTXReader, FASTXTranslator


# Sample files
SAMPLE_FILES = {
    'GZIPPABLE':
        b'>Test\n' + b'ACGATCGATCGATCGAACGATCGTACGTAGCCGTCGATCGACACGA\n' * 30,
    'VALID_FASTQ':
        (b'@Header1\nACGTACGTACGT\n+Header1\nAAAAAAAAAAAA\n'
         b'@Header2\nACGTACGTACGT\n+Header2\nAAAAAAAAAAAA\n'),
    'INVALID_FASTQ':
        (b'@Header1\nACGTACGTACGT\n+\nAAAAAAAAAAAA\n'
         b'@Header2\nACGTACGTAQGT\n+\nAAAAAAAAAAAA\n'),
    'MODIFIABLE_FASTQ':
        (b'@Header1\nACGTACGTACGT\n+\nAAAAAAAAAAAA\n'
         b'@Header2\nACGTACGTAXGT\n+\nAAAAAAAAAAAA\n'),
    'TABBED_FASTQ':
        (b'@Header1\tSecret Val\nACGTACGTACGT\n+\nAAAAAAAAAAAA\n'
         b'@Header2\nACGTACGTANGT\n+\nAAAAAAAAAAAA\n'),
}


def test_nucl_iterator():
    fakefile = BytesIO(b'>test\nACGT\n')
    iterator = iter(FASTXNuclIterator(fakefile))
    val = next(iterator)
    assert val == b'>test\nACGT\n'
    with pytest.raises(StopIteration):
        next(iterator)


def test_paired_validator():
    # test a single file
    fakefile = BytesIO(b'>test\nACGT\n')
    outfile = FASTXTranslator(fakefile, recompress=False)
    assert outfile.read() == b'>test\nACGT\n'

    # test a single file without an ending newline
    fakefile = BytesIO(b'>test\nACGT')
    outfile = FASTXTranslator(fakefile, recompress=False)
    assert outfile.read() == b'>test\nACGT\n'

    # test paired files
    fakefile = BytesIO(b'>test\nACGT\n')
    fakefile2 = BytesIO(b'>test2\nTGCA\n')
    outfile = FASTXTranslator(fakefile, pair=fakefile2, recompress=False)
    assert outfile.read() == b'>test\nACGT\n>test2\nTGCA\n'

    # test compression works
    fakefile = BytesIO(b'>test\nACGT\n')
    outfile = FASTXTranslator(fakefile)
    outdata = outfile.read()

    # there's a 4-byte timestamp in the middle of the gziped data so we check the start and end
    assert outdata.startswith(b'\x1f\x8b\x08\x00')
    assert outdata.endswith(b'\x02\xff\xb3+I-.\xe1rtv\x0f\xe1\x02\x00\xf3\x1dK\xc4\x0b\x00\x00\x00')


def test_file_size_requirement(runner):
    # File must be >= 70 bytes
    with runner.isolated_filesystem():
        with open('myfasta.fa', mode='w') as f:
            f.write('>Too Small\nACGT\n')
        with pytest.raises(ValidationError):
            FASTXTranslator(open('myfasta.fa'))


def test_gzip_filename_validation(runner):
    with runner.isolated_filesystem():
        with gzip.open('myfasta.fa.gz.extra', mode='w') as f:
            f.write(SAMPLE_FILES['GZIPPABLE'])
        with pytest.raises(ValidationError):
            FASTXTranslator(open('myfasta.fa.gz.extra', mode='rb'))


def test_mislabeled_gzip_filename(runner):
    with runner.isolated_filesystem():
        with open('myfasta.fa.gz', mode='w') as f:
            f.write(str(SAMPLE_FILES['GZIPPABLE']))
        with pytest.raises(ValidationError):
            FASTXTranslator(open('myfasta.fa.gz', mode='rb'))


@pytest.mark.skipif(sys.version_info < (3, 3),
                    reason="bz2.open interface requires python3.3")
def test_bzip_filename_validation(runner):
    with runner.isolated_filesystem():
        with bz2.BZ2File('myfasta.fa.bz2.extra', mode='w') as f:
            f.write(SAMPLE_FILES['GZIPPABLE'])
        with pytest.raises(ValidationError):
            FASTXTranslator(open('myfasta.fa.bz2.extra', mode='rb'))


@pytest.mark.skipif(sys.version_info < (3, 3),
                    reason="bz2.open interface requires python3.3")
def test_bzip_file(runner):
    with runner.isolated_filesystem():
        with bz2.open('myfasta.fa.bz2', mode='w') as f:
            f.write(SAMPLE_FILES['GZIPPABLE'])
        translator = FASTXTranslator(open('myfasta.fa.bz2', mode='rb'), recompress=False)
        assert translator.read() == SAMPLE_FILES['GZIPPABLE']


def test_translator_to_reader(runner):
    with runner.isolated_filesystem():
        with gzip.open('myfasta.fa.gz', mode='w') as f:
            f.write(SAMPLE_FILES['GZIPPABLE'])

        translator = FASTXTranslator(open('myfasta.fa.gz', mode='rb'))
        translator.validate()
        with pytest.raises(AssertionError):
            translator.close()
        reader = FASTXReader(translator.reads.file_obj.fileobj)
        assert reader.reads.tell() == 0
        assert reader.readall() == open('myfasta.fa.gz', mode='rb').read()  # Compressed!
        reader.seek(0)
        assert reader.reads.tell() == 0
        reader.close()


@pytest.mark.parametrize('file_id,filename,validates,allow_iupac,modified', [
    ('VALID_FASTQ', 'my.fq', True, False, False),
    ('INVALID_FASTQ', 'my.fq', False, False, False),
    ('MODIFIABLE_FASTQ', 'my.fq', True, True, True),
    ('TABBED_FASTQ', 'my.fq', True, False, True),
    ('VALID_FASTQ', 'my.fastz', False, False, False),  # Bad name
    ('VALID_FASTQ', 'my.fq.gz', False, False, False),   # Bad extension, not gzipped
])
def test_fastq_validation(runner, filename, file_id, validates, allow_iupac, modified):
    warnings.filterwarnings('ignore', category=ValidationWarning)
    with runner.isolated_filesystem():
        with open(filename, mode='wb') as f:
            f.write(SAMPLE_FILES[file_id])

        if validates:
            translator = FASTXTranslator(open(filename, mode='rb'),
                                         allow_iupac=allow_iupac, recompress=False)
            content = SAMPLE_FILES[file_id].replace(b'X', b'N').replace(b'\t', b'|')
            assert translator.read() == content
            assert translator.modified == modified
        else:
            with pytest.raises(ValidationError):
                translator = FASTXTranslator(open(filename, mode='rb'),
                                             allow_iupac=allow_iupac, recompress=False)
                translator.read()


def test_gzip_correctness_large_file(runner):
    """Test that a large random file is properly written out and read back in
    """
    content = '>test_{}\n'
    file_content = ''
    with runner.isolated_filesystem():
        N_READS = 5000
        with open('myfasta.fa', mode='w') as f:
            for ix in range(N_READS):
                read = content.format(ix)
                read += ''.join(random.choice('ACGT') for _ in range(100)) + '\n'
                f.write(read)
                file_content += read

        outfile = FASTXTranslator(open('myfasta.fa', 'rb'), recompress=True)
        with open('myfasta.fa.gz', mode='wb') as f:
            f.write(outfile.read())

        # Ensure we're actually testing the buffer flushing
        assert len(file_content) > 2 * outfile.checked_buffer.MAX_READS_BUFFER_SIZE

        # Check that content matches
        translated_content = gzip.open('myfasta.fa.gz', 'rb').read()
        assert translated_content == file_content.encode()


@pytest.mark.parametrize('n_newlines', [0, 1, 2, 3])
def test_validator_newlines(runner, n_newlines):
    # Test that we properly parse files with 0, 1, or more than 1 newline at the end
    # including our n bytes remaining calculations
    # (Use runner for isolated filesystem only here)
    content = (b'>test\nACGTACGTAGCTGACACGTACGTAGCTGACACGTACGTAGCTGACACGTACGTAGCTGAC'
               b'ACGTACGTAGCTGACACGTACGTAGCTGACACGTACGTAGCTGACACGTACGTAGCTGAC')
    content += (b'\n' * n_newlines)
    fakefile = BytesIO(content)
    outfile = FASTXTranslator(fakefile, recompress=False)
    if n_newlines == 0:
        assert content[-1] == 'C'.encode()[0]
    else:
        assert content[(-1 * n_newlines):] == b'\n' * n_newlines
    assert outfile.read() == content.decode().rstrip().encode() + b'\n'
    assert outfile.reads.bytes_left == 0

    with runner.isolated_filesystem():
        with open('myfasta.fa', mode='wb') as f:
            f.write(content)
        outfile2 = FASTXTranslator(open('myfasta.fa', 'rb'), recompress=False)
        assert outfile2.read() == content.decode().rstrip().encode() + b'\n'
        assert outfile2.reads.bytes_left == 0


def test_wrapper_speed():
    # to check serialization speed run just this test with:
    # %timeit !py.test tests/test_upload.py::test_wrapper_speed

    # most of the current bottleneck appears to be in the gzip.write step
    long_seq = b'ACGT' * 20
    data = b'\n'.join(b'>header_' + str(i).encode() + b'\n' + long_seq + b'\n'
                      for i in range(200))
    wrapper = FASTXTranslator(BytesIO(data))
    assert len(wrapper.read()) < len(data)
