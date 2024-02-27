import mock
import pytest
from contextlib import contextmanager

pytest.importorskip("pandas")  # noqa

import numpy as np
import pandas as pd

from onecodex.models.collection import SampleCollection
from onecodex.exceptions import OneCodexException, StatsException, StatsWarning
from onecodex.lib.enums import AlphaDiversityStatsTest
from onecodex.stats import AlphaDiversityStatsResults


@contextmanager
def does_not_raise():
    yield


@pytest.mark.parametrize("alpha", [-1, 0, 1, 1.1, 2])
def test_alpha_diversity_stats_invalid_alpha(samples, alpha):
    with pytest.raises(StatsException, match="`alpha` must be between 0 and 1"):
        samples.alpha_diversity_stats(group_by="wheat", alpha=alpha)


@pytest.mark.parametrize("fields", [1, 4.2, False, b"", (42,), [-1.0], ("wheat", 43)])
def test_alpha_diversity_stats_invalid_metadata_fields(samples, fields):
    with pytest.raises(OneCodexException, match="type str"):
        samples.alpha_diversity_stats(group_by=fields)

    with pytest.raises(OneCodexException, match="type str"):
        samples.alpha_diversity_stats(group_by="wheat", paired_by=fields)


def test_alpha_diversity_stats_missing_alpha_diversity_values(samples):
    samples.metadata["col"] = ["a", "b", "a"]

    with pytest.warns(StatsWarning), pytest.raises(
        StatsException, match="`group_by` must have at least two groups.*found 0"
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
    "group_by,num_groups", [("col1", 1), (("col1", "col2"), 0), (["col1", "col2"], 0)]
)
def test_alpha_diversity_stats_missing_group_by_data(samples, group_by, num_groups):
    samples.metadata["col1"] = ["a", None, "a"]
    samples.metadata["col2"] = ["1", "42", np.nan]

    with pytest.warns(StatsWarning), pytest.raises(
        StatsException, match=f"`group_by` must have at least two groups.*found {num_groups}"
    ):
        samples.alpha_diversity_stats(group_by=group_by)


@pytest.mark.parametrize(
    "paired_by,num_groups", [("col1", 1), (("col1", "col2"), 0), (["col1", "col2"], 0)]
)
def test_alpha_diversity_stats_missing_paired_by_data(samples, paired_by, num_groups):
    samples.metadata["group"] = "g1"
    samples.metadata["col1"] = ["a", None, "a"]
    samples.metadata["col2"] = ["1", "42", np.nan]

    with pytest.warns(StatsWarning), pytest.raises(
        StatsException, match=f"`group_by` must have at least two groups.*found {num_groups}"
    ):
        samples.alpha_diversity_stats(group_by="group", paired_by=paired_by)


@pytest.mark.parametrize(
    "values,expectation,num_groups",
    [
        ("one_value", does_not_raise(), 1),
        (["all", "different", "values"], pytest.warns(StatsWarning), 0),
    ],
)
def test_alpha_diversity_stats_not_enough_groups(samples, values, expectation, num_groups):
    samples.metadata["col"] = values

    with expectation, pytest.raises(
        StatsException, match=f"`group_by` must have at least two groups.*found {num_groups}"
    ):
        samples.alpha_diversity_stats(group_by="col")


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
            "`group_by` must have exactly two groups",
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
        group_by_variable="group",
        groups={"g1", "g2"},
        paired_by_variable="pair",
    )


def test_mannwhitneyu_wrong_number_of_groups(samples):
    df = pd.DataFrame({"group": ["g1", "g2", "g3", "g1", "g2"], "shannon": [1, 2, 3, 4, 5]})

    with pytest.raises(StatsException, match="`group_by` must have exactly two groups"):
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
    assert results.group_by_variable == "group"
    assert results.groups == {"g1", "g2", "g3"}
    assert results.posthoc_df is None

    results = samples._kruskal(df, "shannon", "group", 0.2)
    labels = pd.Index(["g1", "g2", "g3"])
    assert (results.posthoc_df.index == labels).all()
    assert (results.posthoc_df.columns == labels).all()
