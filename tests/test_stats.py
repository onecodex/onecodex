import mock
import pytest

pytest.importorskip("pandas")  # noqa

import numpy as np
import pandas as pd
import skbio

from onecodex.models.collection import SampleCollection
from onecodex.exceptions import OneCodexException, StatsException, StatsWarning
from onecodex.lib.enums import AlphaDiversityStatsTest, BetaDiversityStatsTest, AlphaDiversityMetric
from onecodex.stats import AlphaDiversityStatsResults, PosthocResults


@pytest.mark.parametrize("method", ["alpha_diversity_stats", "beta_diversity_stats"])
@pytest.mark.parametrize("alpha", [-1, 0, 1, 1.1, 2])
def test_invalid_alpha(samples, method, alpha):
    with pytest.raises(StatsException, match="`alpha` must be between 0 and 1"):
        getattr(samples, method)(group_by="wheat", alpha=alpha)


@pytest.mark.parametrize("method", ["alpha_diversity_stats", "beta_diversity_stats"])
@pytest.mark.parametrize("fields", [1, 4.2, False, b"", (42,), [-1.0], ("wheat", 43)])
def test_invalid_metadata_fields(samples, method, fields):
    with pytest.raises(OneCodexException, match="type str"):
        getattr(samples, method)(group_by=fields)


@pytest.mark.parametrize("method", ["alpha_diversity_stats", "beta_diversity_stats"])
@pytest.mark.parametrize(
    "group_by,num_groups", [("col1", 1), (("col1", "col2"), 0), (["col1", "col2"], 0)]
)
def test_missing_group_by_data(samples, method, group_by, num_groups):
    samples.metadata["col1"] = ["a", None, "a"]
    samples.metadata["col2"] = ["1", "42", np.nan]

    with (
        pytest.warns(StatsWarning),
        pytest.raises(
            StatsException, match=f"`group_by` must have at least 2 groups.*found {num_groups}"
        ),
    ):
        getattr(samples, method)(group_by=group_by)


@pytest.mark.parametrize(
    "method,kwargs",
    [
        ("alpha_diversity_stats", {"metric": AlphaDiversityMetric.ObservedTaxa}),
        ("beta_diversity_stats", {}),
    ],
)
def test_samples_missing_abundances(samples, method, kwargs):
    samples.metadata["col"] = ["a", "b", "a"]
    samples._results[:] = np.nan
    assert len(samples._classification_ids_without_abundances) == 3

    with (
        pytest.warns(StatsWarning, match="3 samples.*missing abundance"),
        pytest.raises(StatsException, match="`group_by` must have at least 2 groups.*found 0"),
    ):
        getattr(samples, method)(group_by="col", **kwargs)


@pytest.mark.parametrize("method", ["alpha_diversity_stats", "beta_diversity_stats"])
@pytest.mark.parametrize(
    "values,num_groups", [("one_value", 1), (["all", "different", "values"], 0)]
)
@pytest.mark.filterwarnings("ignore:.*size < 2.*")
def test_not_enough_groups(samples, method, values, num_groups):
    samples.metadata["col"] = values

    with pytest.raises(
        StatsException, match=f"`group_by` must have at least 2 groups.*found {num_groups}"
    ):
        getattr(samples, method)(group_by="col")


@pytest.mark.parametrize("paired_by", [1, 4.2, False, b"", (42,), [-1.0], ("wheat", 43)])
def test_alpha_diversity_stats_invalid_paired_by(samples, paired_by):
    with pytest.raises(OneCodexException, match="type str"):
        samples.alpha_diversity_stats(group_by="wheat", paired_by=paired_by)


def test_alpha_diversity_stats_missing_alpha_diversity_values(samples):
    samples.metadata["col"] = ["a", "b", "a"]

    with (
        pytest.warns(StatsWarning),
        pytest.raises(StatsException, match="`group_by` must have at least 2 groups.*found 0"),
    ):
        with mock.patch.object(
            SampleCollection,
            "alpha_diversity",
            return_value=pd.DataFrame(
                {"shannon": [np.nan, np.nan, np.nan]}, index=samples.to_df().index
            ),
        ):
            samples.alpha_diversity_stats(group_by="col")


@pytest.mark.parametrize(
    "paired_by,num_groups", [("col1", 1), (("col1", "col2"), 0), (["col1", "col2"], 0)]
)
def test_alpha_diversity_stats_missing_paired_by_data(samples, paired_by, num_groups):
    samples.metadata["group"] = "g1"
    samples.metadata["col1"] = ["a", None, "a"]
    samples.metadata["col2"] = ["1", "42", np.nan]

    with (
        pytest.warns(StatsWarning),
        pytest.raises(
            StatsException, match=f"`group_by` must have at least 2 groups.*found {num_groups}"
        ),
    ):
        samples.alpha_diversity_stats(group_by="group", paired_by=paired_by)


@pytest.mark.parametrize(
    "groups,paired_by,test",
    [
        (["g1", "g2", "g1", "g2", "g2", "g2"], None, AlphaDiversityStatsTest.Mannwhitneyu),
        (["g1", "g2", "g1", "g2", "g1", "g2"], "pair", AlphaDiversityStatsTest.Wilcoxon),
        (["g1", "g2", "g3", "g1", "g2", "g3"], None, AlphaDiversityStatsTest.Kruskal),
    ],
)
def test_alpha_diversity_stats_auto_test(ocx, api_data, samples, groups, paired_by, test):
    with pytest.warns(UserWarning):
        samples.extend(
            [
                ocx.Samples.get("cc18208d98ad48b3"),
                ocx.Samples.get("5445740666134eee"),
                ocx.Samples.get("0ecac25ec0004fe4"),
            ]
        )
        samples.metadata["group"] = groups
        samples.metadata["pair"] = ["p1", "p1", "p2", "p2", "p3", "p3"]

        results = samples.alpha_diversity_stats(group_by="group", paired_by=paired_by)

        assert results.test == test


@pytest.mark.parametrize(
    "test", [AlphaDiversityStatsTest.Mannwhitneyu, AlphaDiversityStatsTest.Kruskal]
)
def test_alpha_diversity_stats_paired_by_with_wrong_test(ocx, api_data, samples, test):
    samples.extend([ocx.Samples.get("cc18208d98ad48b3")])
    samples.metadata["group"] = ["g1", "g2", "g1", "g2"]
    samples.metadata["pair"] = ["p1", "p1", "p2", "p2"]

    with pytest.raises(StatsException, match="paired_by.*wilcoxon"):
        samples.alpha_diversity_stats(group_by="group", paired_by="pair", test=test)


def test_wilcoxon_missing_paired_by(samples):
    df = pd.DataFrame({"group": ["g1", "g2", "g1", "g2"], "shannon": [1, 2, 3, 4]})

    with pytest.raises(StatsException, match="paired_by.*wilcoxon"):
        samples._wilcoxon(df, "shannon", "group", None)


@pytest.mark.parametrize(
    "group,pair,metric,match",
    [
        # Wrong number of groups
        (
            ["g1", "g2", "g3"],
            ["p1", "p1", "p2"],
            [1, 2, 3],
            "`group_by` must have exactly 2 groups",
        ),
        # Different group sizes
        (
            ["g1", "g2", "g1", "g2", "g2"],
            ["p1", "p1", "p2", "p2", "p2"],
            [1, 2, 3, 4, 5],
            "two groups of the same size",
        ),
        # Duplicate paired by
        (
            ["g1", "g2", "g1", "g2"],
            ["p1", "p1", "p2", "p1"],
            [1, 2, 3, 4],
            "unique values within each group",
        ),
        # Unpaired samples
        (
            ["g1", "g2", "g1", "g2"],
            ["p1", "p1", "p2", "p3"],
            [1, 2, 3, 4],
            "matching paired sample",
        ),
    ],
)
def test_wilcoxon_errors(samples, group, pair, metric, match):
    df = pd.DataFrame({"group": group, "pair": pair, "metric": metric})

    with pytest.raises(StatsException, match=match):
        samples._wilcoxon(df, "metric", "group", "pair")


def test_wilcoxon(samples):
    df = pd.DataFrame(
        {
            "group": ["g1", "g2", "g2", "g1", "g1", "g2"],
            "pair": ["p1", "p2", "p1", "p2", "p3", "p3"],
            "metric": [42.0, 0.0, 45.7, 1.2, 2.3, 4.5],
        }
    )

    results = samples._wilcoxon(df, "metric", "group", "pair")

    # Results verified manually with `scipy.stats.wilcoxon`
    assert results == AlphaDiversityStatsResults(
        test=AlphaDiversityStatsTest.Wilcoxon,
        statistic=1.0,
        pvalue=0.5,
        sample_size=6,
        group_by_variable="group",
        groups={"g1", "g2"},
        paired_by_variable="pair",
    )


def test_mannwhitneyu_wrong_number_of_groups(samples):
    df = pd.DataFrame({"group": ["g1", "g2", "g3", "g1", "g2"], "shannon": [1, 2, 3, 4, 5]})

    with pytest.raises(StatsException, match="`group_by` must have exactly 2 groups"):
        samples._mannwhitneyu(df, "shannon", "group")


def test_mannwhitneyu(samples):
    df = pd.DataFrame(
        {"group": ["g1", "g2", "g1", "g2", "g2"], "shannon": [42.0, 0.0, 45.7, 1.2, 2.3]}
    )

    results = samples._mannwhitneyu(df, "shannon", "group")

    # Results verified manually with `scipy.stats.mannwhitneyu`
    assert results == AlphaDiversityStatsResults(
        test=AlphaDiversityStatsTest.Mannwhitneyu,
        statistic=6.0,
        pvalue=0.2,
        sample_size=5,
        group_by_variable="group",
        groups={"g1", "g2"},
    )


def test_kruskal(samples):
    df = pd.DataFrame(
        {"group": ["g1", "g2", "g3", "g2", "g1", "g3"], "shannon": [42.0, 0.0, 45.7, 1.2, 2.3, 4.5]}
    )

    results = samples._kruskal(df, "shannon", "group", 0.05)

    # Results verified manually with `scipy.stats.kruskal`
    assert results.test == AlphaDiversityStatsTest.Kruskal
    assert results.statistic.round(3) == 3.714
    assert results.pvalue.round(3) == 0.156
    assert results.sample_size == 6
    assert results.group_by_variable == "group"
    assert results.groups == {"g1", "g2", "g3"}
    assert results.posthoc is None

    results = samples._kruskal(df, "shannon", "group", 0.2)
    labels = pd.Index(["g1", "g2", "g3"])
    assert (results.posthoc.adjusted_pvalues.index == labels).all()
    assert (results.posthoc.adjusted_pvalues.columns == labels).all()


def test_alpha_diversity_stats_results_posthoc_df_property():
    results = AlphaDiversityStatsResults(
        test=AlphaDiversityStatsTest.Mannwhitneyu,
        statistic=6.0,
        pvalue=0.2,
        sample_size=5,
        group_by_variable="group",
        groups={"g1", "g2"},
        posthoc=None,
    )

    with pytest.warns(DeprecationWarning):
        assert results.posthoc_df is None

    pvals = pd.DataFrame([[1.0, 0.2], [0.2, 1.0]], index=["g1", "g2"], columns=["g1", "g2"])
    results = AlphaDiversityStatsResults(
        test=AlphaDiversityStatsTest.Mannwhitneyu,
        statistic=6.0,
        pvalue=0.2,
        sample_size=5,
        group_by_variable="group",
        groups={"g1", "g2"},
        posthoc=PosthocResults(adjusted_pvalues=pvals),
    )

    with pytest.warns(DeprecationWarning):
        assert results.posthoc_df.equals(pvals)


def test_permanova(samples):
    ids = ["id1", "id2", "id3", "id4", "id5", "id6"]
    dm = skbio.DistanceMatrix(
        [
            [0, 1, 20, 42, 2, 3],
            [1, 0, 38, 24, 4, 5],
            [20, 38, 0, 1.5, 10, 20],
            [42, 24, 1.5, 0, 12, 11],
            [2, 4, 10, 12, 0, 1],
            [3, 5, 20, 11, 1, 0],
        ],
        ids,
    )
    df = pd.DataFrame({"group": ["g1", "g1", "g2", "g2", "g3", "g3"]}, index=ids)

    np.random.seed(0)  # deterministic p-values
    results = samples._permanova(dm, df, "group", 0.04, 999)

    # Results verified manually with `skbio.stats.distance.permanova`
    assert results.test == BetaDiversityStatsTest.Permanova
    assert results.statistic.round(3) == 587.588
    assert results.pvalue.round(3) == 0.048
    assert results.num_permutations == 999
    assert results.sample_size == 6
    assert results.group_by_variable == "group"
    assert results.groups == {"g1", "g2", "g3"}
    assert results.posthoc is None

    np.random.seed(0)
    posthoc = samples._permanova(dm, df, "group", 0.05, 999).posthoc
    labels = pd.Index(["g1", "g2", "g3"])
    for attr in "statistics", "pvalues", "adjusted_pvalues":
        df = getattr(posthoc, attr)
        assert (df.index == labels).all()
        assert (df.columns == labels).all()

        data = df.values
        # Data is symmetric (ignoring the possibly NaN diagonal)
        upper_triangle_indices = np.triu_indices_from(data, k=1)
        assert (data[upper_triangle_indices] == data.T[upper_triangle_indices]).all()

        # ...with correct diagonal and upper triangle values
        if attr == "statistics":
            assert np.isnan(np.diagonal(data)).all()

            statistics = [1286.385, 26.0, 234.385]
            assert (np.round(data[upper_triangle_indices], 3) == statistics).all()
        elif attr == "pvalues":
            assert (np.diagonal(data) == 1.0).all()

            pvalues = [0.351, 0.315, 0.352]
            assert (np.round(data[upper_triangle_indices], 3) == pvalues).all()
        elif attr == "adjusted_pvalues":
            assert (np.diagonal(data) == 1.0).all()

            adjusted_pvalues = [0.352, 0.352, 0.352]
            assert (np.round(data[upper_triangle_indices], 3) == adjusted_pvalues).all()
