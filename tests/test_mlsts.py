from onecodex.models import Mlsts


def test_query_for_mlsts(ocx, api_data):
    sample_id = "002f55a11aae4b8f"
    mlst = ocx.Mlsts.where(sample=sample_id)
    assert isinstance(mlst[0], Mlsts)

    profiles = ocx.Mlsts.all()
    assert len(profiles) == 1


def test_mlsts_results(ocx, api_data):
    mlst_id = "cff9548162fa4c9b"
    mlst = Mlsts.get(mlst_id)

    assert isinstance(mlst, Mlsts)

    result = mlst.results()

    assert "mlstResults" in result
    assert result["mlstResults"]["summary"]["mlst"] == 42
    assert result["mlstResults"]["results"][0]["st"] == 42
