import pytest

pytest.importorskip("pandas")  # noqa
import pandas as pd
import skbio

from onecodex.exceptions import OneCodexException


@pytest.mark.parametrize(
    "metric,value,kwargs",
    [
        ("simpson", [0.9194601552055255, 0.8918602435339009, 0.7805781921054538], {}),
        ("chao1", [157.0, 128.0, 97.0], {}),
    ],
)
def test_alpha_diversity(ocx, api_data, metric, value, kwargs):
    samples = ocx.Samples.where(project="4b53797444f846c4")
    divs = samples.alpha_diversity(metric=metric, **kwargs)
    assert isinstance(divs, pd.DataFrame)
    assert divs[metric].tolist() == value


def test_alpha_diversity_exceptions(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    # must be a metric that exists
    with pytest.raises(OneCodexException) as e:
        samples.alpha_diversity("does_not_exist")
    assert "metric must be one of" in str(e.value)


@pytest.mark.parametrize(
    "metric,value,kwargs",
    [
        ("braycurtis", [0.886539, 0.854311, 0.909196], {}),
        ("cityblock", [1.773079, 1.708622, 1.818392], {}),
        ("jaccard", [0.755459, 0.71066, 0.770492], {}),
        ("unweighted_unifrac", [0.618799, 0.562130, 0.615152], {}),
        ("weighted_unifrac", [0.507442, 0.405466, 0.612872], {}),
    ],
)
def test_beta_diversity(ocx, api_data, metric, value, kwargs):
    samples = ocx.Samples.where(project="4b53797444f846c4")
    dm = samples.beta_diversity(metric=metric, **kwargs)
    assert isinstance(dm, skbio.stats.distance._base.DistanceMatrix)
    assert dm.condensed_form().round(6).tolist() == value


def test_beta_diversity_exceptions(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    # must be a metric that exists
    with pytest.raises(OneCodexException) as e:
        samples.beta_diversity("does_not_exist")
    assert "metric must be one of" in str(e.value)


@pytest.mark.parametrize(
    "weighted,value",
    [(False, [0.618799, 0.562130, 0.615152]), (True, [0.507442, 0.405466, 0.612872])],
)
def test_unifrac(ocx, api_data, value, weighted):
    samples = ocx.Samples.where(project="4b53797444f846c4")
    dm = samples.unifrac(weighted=weighted)
    assert isinstance(dm, skbio.stats.distance._base.DistanceMatrix)
    assert dm.condensed_form().round(6).tolist() == value
