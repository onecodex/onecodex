import pytest

pytest.importorskip("pandas")  # noqa
import pandas as pd
import skbio

from onecodex.exceptions import OneCodexException


@pytest.mark.parametrize(
    "metric,value,kwargs",
    [
        ("simpson", [0.9232922257199748, 0.8930761430647977, 0.7865654458730155], {}),
        ("chao1", [164.0, 134.0, 103.0], {}),
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
        ("braycurtis", [0.886014, 0.84694, 0.905716], {}),
        ("cityblock", [1.772028, 1.693879, 1.811432], {}),
        ("jaccard", [0.742616, 0.697561, 0.752632], {}),
        ("unweighted_unifrac", [0.6, 0.547486, 0.591304], {}),
        ("weighted_unifrac", [0.503168, 0.403155, 0.605155], {}),
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
    "weighted,value", [(False, [0.6, 0.547486, 0.591304]), (True, [0.503168, 0.403155, 0.605155])],
)
def test_unifrac(ocx, api_data, value, weighted):
    samples = ocx.Samples.where(project="4b53797444f846c4")
    dm = samples.unifrac(weighted=weighted)
    assert isinstance(dm, skbio.stats.distance._base.DistanceMatrix)
    assert dm.condensed_form().round(6).tolist() == value
