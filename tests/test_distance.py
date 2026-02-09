import mock
import pytest

pytest.importorskip("pandas")  # noqa
import pandas as pd
import skbio

from onecodex.exceptions import OneCodexException
from onecodex.lib.enums import BetaDiversityMetric


@pytest.mark.parametrize(
    "metric,value",
    [
        ("simpson", [0.9232922257199748, 0.8930761430647977, 0.7865654458730156]),
        ("observed_taxa", [164.0, 134.0, 103.0]),
    ],
)
def test_alpha_diversity(samples, metric, value):
    divs = samples.alpha_diversity(diversity_metric=metric, metric="abundance_w_children")
    assert isinstance(divs, pd.DataFrame)
    assert divs[metric].tolist() == value


def test_alpha_diversity_exceptions(samples):
    # must be a metric that exists
    with pytest.raises(OneCodexException) as e:
        samples.alpha_diversity(diversity_metric="does_not_exist")
    assert "metric must be one of" in str(e.value)


def test_alpha_diversity_warnings(samples):
    with pytest.warns(DeprecationWarning, match="`Chao1` is deprecated"):
        samples.alpha_diversity(diversity_metric="chao1")


@pytest.mark.parametrize(
    "metric,value,kwargs",
    [
        ("braycurtis", [0.886014, 0.84694, 0.905716], {}),
        ("cityblock", [1.772028, 1.693879, 1.811432], {}),
        ("manhattan", [1.772028, 1.693879, 1.811432], {}),
        ("jaccard", [0.742616, 0.697561, 0.752632], {}),
        ("unweighted_unifrac", [0.6, 0.547486, 0.591304], {}),
        ("weighted_unifrac", [0.503168, 0.403155, 0.605155], {}),
        ("aitchison", [59.800677, 51.913911, 52.693709], {}),
    ],
)
def test_beta_diversity(samples, metric, value, kwargs):
    assert samples.automatic_rank(metric="auto") == "species"
    dm = samples.beta_diversity(rank="species", diversity_metric=metric, **kwargs)
    assert isinstance(dm, skbio.stats.distance._base.DistanceMatrix)
    assert dm.condensed_form().round(6).tolist() == value


def test_beta_diversity_braycurtis_nans(samples):
    mock_df = samples.to_df()
    mock_df.loc[:, :] = 0

    with mock.patch.object(samples, "to_df", return_value=mock_df):
        dm = samples.beta_diversity(diversity_metric=BetaDiversityMetric.BrayCurtis)

    assert isinstance(dm, skbio.stats.distance._base.DistanceMatrix)
    assert (dm.condensed_form() == 0.0).all()


def test_beta_diversity_exceptions(samples):
    # must be a metric that exists
    with pytest.raises(OneCodexException) as e:
        samples.beta_diversity(diversity_metric="does_not_exist")
    assert "metric must be one of" in str(e.value)


@pytest.mark.parametrize(
    "weighted,value",
    [(False, [0.6, 0.547486, 0.591304]), (True, [0.503168, 0.403155, 0.605155])],
)
def test_unifrac(samples, value, weighted):
    dm = samples.unifrac(weighted=weighted)
    assert isinstance(dm, skbio.stats.distance._base.DistanceMatrix)
    assert dm.condensed_form().round(6).tolist() == value


# This tests that Unifrac calculations don't break when `root` has
# more than one child by mocking a new node that is a direct child
# of `root`. There's a bug in scikit-bio where it requires that
# the root has only one child, which isn't true in our taxonomy.
# See onecodex/distances.py for more details.
def test_unifrac_tree(samples):
    samples[0].primary_classification.results()  # seed _cached_results
    samples[0].primary_classification._cached_result["table"].append(
        {
            "abundance": None,
            "name": "fake node",
            "parent_tax_id": "1",
            "rank": "species",
            "readcount": 100000,
            "readcount_w_children": 100000,
            "tax_id": "1000000000",
        }
    )

    df = samples.to_df(metric="readcount_w_children")
    assert df.shape[1] == 1081  # make sure our insert worked

    dm = samples.unifrac()
    assert isinstance(dm, skbio.stats.distance._base.DistanceMatrix)
