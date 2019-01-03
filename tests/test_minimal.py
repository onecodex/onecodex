import pytest


def test_minimal():
    with pytest.raises(ImportError):
        import altair  # noqa
