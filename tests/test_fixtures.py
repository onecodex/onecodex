from onecodex import Api


def test_api_fixture(ocx):
    assert isinstance(ocx, Api)
