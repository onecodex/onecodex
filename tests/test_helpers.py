from onecodex.helpers import is_insecure_platform


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
        assert is_insecure_platform() is True
