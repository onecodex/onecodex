import pytest


def test_minimal():
    with pytest.raises(ImportError):
        import matplotlib  # noqa
