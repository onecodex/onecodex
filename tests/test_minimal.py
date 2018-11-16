import pytest

import onecodex
from onecodex import Api


def test_minimal():
    with pytest.raises(ImportError):
        import matplotlib  # noqa


def test_modulealias():
    # ensure viz and distance modules aren't imported if matplotlib isn't around
    ocx = Api()

    with pytest.raises(AttributeError):
        isinstance(ocx.viz, onecodex.utils.ModuleAlias)
    with pytest.raises(AttributeError):
        isinstance(ocx.distance, onecodex.utils.ModuleAlias)
