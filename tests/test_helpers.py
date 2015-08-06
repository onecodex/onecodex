from onecodex.helpers import check_for_allowed_file, is_insecure_platform


def test_valid_extensions():
    files = [
        "my.fasta",
        "my.fastq",
        "my.fasta.gz",
        "my.fastq.gz",
        "my.fa",
        "my.fq",
        "my.fa.gz",
        "my.fq.gz",
        "my.really.long.file.fasta",
        "odd-dashes-in-file.fa"
    ]
    for f in files:
        assert check_for_allowed_file(f)


def test_warnings():
    import sys
    v = sys.version_info
    if v.major == 2 and v.minor == 7 and v.micro >= 9:
        assert is_insecure_platform() is False

    try:
        import OpenSSL  # noqa
        import ndg  # noqa
        import pyasn1  # noqa
        assert is_insecure_platform() is False
    except ImportError:
        import ssl  # noqa
        if v.major == 2 and v.minor == 7 and v.micro >= 9:
            assert is_insecure_platform() is False
        else:
            assert is_insecure_platform() is True
