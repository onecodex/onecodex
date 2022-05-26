from onecodex import Api


def test_api_fixture(ocx):
    assert isinstance(ocx, Api)


def test_experimental_api_fixture(ocx_experimental):
    assert isinstance(ocx_experimental, Api)
    assert "v1_experimental" in ocx_experimental._schema_path
