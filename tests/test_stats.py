import mock
import pytest

pytest.importorskip("pandas")  # noqa

import numpy as np
import pandas as pd
import skbio

from onecodex.exceptions import OneCodexException, StatsException, StatsWarning
from onecodex.lib.enums import (
    AdjustmentMethod,
    AlphaDiversityMetric,
    AlphaDiversityStatsTest,
    BetaDiversityStatsTest,
    Metric,
    PosthocStatsTest,
)
from onecodex.models.collection import SampleCollection
from onecodex.stats import (
    AlphaDiversityStatsResults,
    AncombcResults,
    BetaDiversityStatsResults,
)

ANCOMBC_MAIN_RESULTS_COLUMNS = {"Log2(FC)", "SE", "W", "pvalue", "qvalue", "Signif"}
ANCOMBC_MAIN_RESULTS_INDEX_NAMES = ["Taxon", "Covariate"]
ANCOMBC_GLOBAL_RESULTS_INDEX_NAME = "Taxon"
ANCOMBC_GLOBAL_RESULTS_COLUMNS = {"W", "pvalue", "qvalue", "Signif"}


@pytest.mark.parametrize("method", ["alpha_diversity_stats", "beta_diversity_stats", "_ancombc"])
@pytest.mark.parametrize("alpha", [-1, 0, 1, 1.1, 2])
def test_invalid_alpha(samples, method, alpha):
    with pytest.raises(StatsException, match="`alpha` must be between 0 and 1"):
        getattr(samples, method)(group_by="wheat", alpha=alpha)


@pytest.mark.parametrize("method", ["alpha_diversity_stats", "beta_diversity_stats", "_ancombc"])
@pytest.mark.parametrize("fields", [1, 4.2, False, b"", (42,), [-1.0], ("wheat", 43)])
def test_invalid_metadata_fields(samples, method, fields):
    with pytest.raises(OneCodexException, match="type str"):
        getattr(samples, method)(group_by=fields)


@pytest.mark.parametrize("method", ["alpha_diversity_stats", "beta_diversity_stats", "_ancombc"])
@pytest.mark.parametrize(
    "group_by,num_groups", [("col1", 1), (("col1", "col2"), 0), (["col1", "col2"], 0)]
)
def test_missing_group_by_data(samples, method, group_by, num_groups):
    for sample, col1, col2 in zip(samples, ["a", None, "a"], ["1", "42", np.nan]):
        sample.metadata.custom["col1"] = col1
        sample.metadata.custom["col2"] = col2

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
        (
            "alpha_diversity_stats",
            {
                "diversity_metric": AlphaDiversityMetric.ObservedTaxa,
                "metric": Metric.AbundanceWChildren,
            },
        ),
        ("beta_diversity_stats", {"metric": Metric.AbundanceWChildren}),
        ("_ancombc", {"metric": Metric.AbundanceWChildren}),
    ],
)
def test_samples_missing_abundances(samples_without_abundances, method, kwargs):
    for sample, col in zip(samples_without_abundances, ["a", "b", "a"]):
        sample.metadata.custom["col"] = col

    with (
        pytest.warns(StatsWarning, match="3 samples.*missing abundance"),
        pytest.raises(StatsException, match="`group_by` must have at least 2 groups.*found 0"),
    ):
        getattr(samples_without_abundances, method)(group_by="col", **kwargs)


@pytest.mark.parametrize(
    "method,kwargs",
    [
        (
            "alpha_diversity_stats",
            {
                "diversity_metric": AlphaDiversityMetric.ObservedTaxa,
                "metric": Metric.ReadcountWChildren,
            },
        ),
        ("beta_diversity_stats", {"metric": Metric.ReadcountWChildren}),
        ("_ancombc", {"metric": Metric.ReadcountWChildren}),
    ],
)
def test_samples_missing_abundance_readcount_metric(samples_without_abundances, method, kwargs):
    """If using a non abundance-based metric; no warning about missing abundances should be issued"""
    for sample, col in zip(samples_without_abundances, ["a", "b", "a"]):
        sample.metadata.custom["col"] = col

    with (
        pytest.warns(StatsWarning, match="1 sample was excluded*"),
        pytest.raises(StatsException, match="`group_by` must have at least 2 groups.*found 1"),
    ):
        getattr(samples_without_abundances, method)(group_by="col", **kwargs)


@pytest.mark.parametrize("method", ["alpha_diversity_stats", "beta_diversity_stats", "_ancombc"])
@pytest.mark.parametrize(
    "values,num_groups",
    [(["one_value", "one_value", "one_value"], 1), (["all", "different", "values"], 0)],
)
@pytest.mark.filterwarnings("ignore:.*size < 2.*")
def test_not_enough_groups(samples, method, values, num_groups):
    for sample, val in zip(samples, values):
        sample.metadata.custom["col"] = val

    with pytest.raises(
        StatsException, match=f"`group_by` must have at least 2 groups.*found {num_groups}"
    ):
        getattr(samples, method)(group_by="col", rank="genus")


@pytest.mark.parametrize("method", ["alpha_diversity_stats", "beta_diversity_stats", "_ancombc"])
@pytest.mark.parametrize("require_classification_version_match", [True, False])
@pytest.mark.filterwarnings("ignore:.*multiple analysis types.*:UserWarning")
def test_require_classification_version_match(
    ocx, api_data, samples, method, require_classification_version_match
):
    # This sample has a different job ID than the others in `samples`
    samples.append(ocx.Samples.get("0ecac25ec0004fe4"))

    # We need to set the metadata in order to reach the error we're trying to trigger
    for sample, val in zip(samples, ["a", "b", "a", "b"]):
        sample.metadata.custom["col"] = val

    kwargs = dict(
        group_by="col",
        rank="genus",
        metric="readcount_w_children",
        require_classification_version_match=require_classification_version_match,
    )

    if require_classification_version_match:
        with pytest.raises(StatsException, match="different database versions"):
            getattr(samples, method)(**kwargs)
    else:
        getattr(samples, method)(**kwargs)


@pytest.mark.parametrize("paired_by", [1, 4.2, False, b"", (42,), [-1.0], ("wheat", 43)])
def test_alpha_diversity_stats_invalid_paired_by(samples, paired_by):
    with pytest.raises(OneCodexException, match="type str"):
        samples.alpha_diversity_stats(group_by="wheat", paired_by=paired_by)


def test_alpha_diversity_stats_missing_alpha_diversity_values(samples):
    for sample, val in zip(samples, ["a", "b", "a"]):
        sample.metadata.custom["col"] = val

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
    for sample, group, col1, col2 in zip(
        samples, ["g1", "g1", "g1"], ["a", None, "a"], ["1", "42", np.nan]
    ):
        sample.metadata.custom["group"] = group
        sample.metadata.custom["col1"] = col1
        sample.metadata.custom["col2"] = col2

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
    samples.extend(
        [
            ocx.Samples.get("cc18208d98ad48b3"),
            ocx.Samples.get("5445740666134eee"),
            # abundance_w_children are all zeroes, which causes alpha diversity to get computed as
            # NaN and the sample to be dropped, so we use readcount_w_children
            ocx.Samples.get("0ecac25ec0004fe4"),
        ]
    )

    pairs = ["p1", "p1", "p2", "p2", "p3", "p3"]

    n_updated = 0
    for sample, group, pair in zip(samples, groups, pairs):
        sample.metadata.custom["group"] = group
        sample.metadata.custom["pair"] = pair
        n_updated += 1

    assert n_updated == len(samples)

    assert list(samples.metadata["group"]) == groups
    assert list(samples.metadata["pair"]) == pairs

    with pytest.warns(UserWarning):  # unrelated -> multiple DBs
        results = samples.alpha_diversity_stats(
            metric="readcount_w_children",
            group_by="group",
            paired_by=paired_by,
            rank="genus",
            require_classification_version_match=False,
        )
    assert results.test == test


@pytest.mark.parametrize(
    "test", [AlphaDiversityStatsTest.Mannwhitneyu, AlphaDiversityStatsTest.Kruskal]
)
def test_alpha_diversity_stats_paired_by_with_wrong_test(ocx, api_data, samples, test):
    samples.extend([ocx.Samples.get("cc18208d98ad48b3")])

    for sample, group, pair in zip(samples, ["g1", "g2", "g1", "g2"], ["p1", "p1", "p2", "p2"]):
        sample.metadata.custom["group"] = group
        sample.metadata.custom["pair"] = pair

    with pytest.raises(StatsException, match="paired_by.*wilcoxon"):
        samples.alpha_diversity_stats(group_by="group", paired_by="pair", test=test, rank="genus")


def test_wilcoxon_missing_paired_by(samples):
    df = pd.DataFrame({"group": ["g1", "g2", "g1", "g2"], "shannon": [1, 2, 3, 4]})

    with pytest.raises(StatsException, match="paired_by.*wilcoxon"):
        samples._wilcoxon(df, "shannon", "group", None, 0.05)


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
        samples._wilcoxon(df, "metric", "group", "pair", 0.05)


def test_wilcoxon(samples):
    df = pd.DataFrame(
        {
            "group": ["g1", "g2", "g2", "g1", "g1", "g2"],
            "pair": ["p1", "p2", "p1", "p2", "p3", "p3"],
            "metric": [42.0, 0.0, 45.7, 1.2, 2.3, 4.5],
        }
    )

    results = samples._wilcoxon(df, "metric", "group", "pair", 0.05)

    # Results verified manually with `scipy.stats.wilcoxon`
    assert results == AlphaDiversityStatsResults(
        test=AlphaDiversityStatsTest.Wilcoxon,
        statistic=1.0,
        pvalue=0.5,
        alpha=0.05,
        sample_size=6,
        group_by_variable="group",
        group_sizes={"g1": 3, "g2": 3},
        paired_by_variable="pair",
    )


def test_mannwhitneyu_wrong_number_of_groups(samples):
    df = pd.DataFrame({"group": ["g1", "g2", "g3", "g1", "g2"], "shannon": [1, 2, 3, 4, 5]})

    with pytest.raises(StatsException, match="`group_by` must have exactly 2 groups"):
        samples._mannwhitneyu(df, "shannon", "group", 0.05)


def test_mannwhitneyu(samples):
    df = pd.DataFrame(
        {"group": ["g1", "g2", "g1", "g2", "g2"], "shannon": [42.0, 0.0, 45.7, 1.2, 2.3]}
    )

    results = samples._mannwhitneyu(df, "shannon", "group", 0.05)

    # Results verified manually with `scipy.stats.mannwhitneyu`
    assert results == AlphaDiversityStatsResults(
        test=AlphaDiversityStatsTest.Mannwhitneyu,
        statistic=6.0,
        pvalue=0.2,
        alpha=0.05,
        sample_size=5,
        group_by_variable="group",
        group_sizes={"g1": 2, "g2": 3},
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
    assert results.alpha == 0.05
    assert results.sample_size == 6
    assert results.group_by_variable == "group"
    assert results.group_sizes == {"g1": 2, "g2": 2, "g3": 2}
    assert results.posthoc is None

    results = samples._kruskal(df, "shannon", "group", 0.2)
    assert results.alpha == 0.2
    assert results.posthoc.test == PosthocStatsTest.Dunn
    assert results.posthoc.adjustment_method == AdjustmentMethod.BenjaminiHochberg
    labels = pd.Index(["g1", "g2", "g3"])
    assert (results.posthoc.adjusted_pvalues.index == labels).all()
    assert (results.posthoc.adjusted_pvalues.columns == labels).all()


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
    assert results.alpha == 0.04
    assert results.num_permutations == 999
    assert results.sample_size == 6
    assert results.group_by_variable == "group"
    assert results.group_sizes == {"g1": 2, "g2": 2, "g3": 2}
    assert results.posthoc is None

    np.random.seed(0)
    results_with_posthoc = samples._permanova(dm, df, "group", 0.05, 999)
    assert results_with_posthoc.alpha == 0.05
    posthoc = results_with_posthoc.posthoc
    assert posthoc.test == PosthocStatsTest.PairwisePermanova
    assert posthoc.adjustment_method == AdjustmentMethod.BenjaminiHochberg
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


@pytest.mark.parametrize(
    "results,expected_groups",
    [
        (
            AlphaDiversityStatsResults(
                test=AlphaDiversityStatsTest.Mannwhitneyu,
                statistic=6.0,
                pvalue=0.2,
                alpha=0.05,
                sample_size=5,
                group_by_variable="group",
                group_sizes={"g1": 2, "g2": 3},
            ),
            {"g1", "g2"},
        ),
        (
            BetaDiversityStatsResults(
                test=BetaDiversityStatsTest.Permanova,
                statistic=1.0,
                pvalue=0.05,
                alpha=0.05,
                num_permutations=999,
                sample_size=6,
                group_by_variable="group",
                group_sizes={"g1": 2, "g2": 2, "g3": 2},
            ),
            {"g1", "g2", "g3"},
        ),
    ],
)
def test_stats_results_groups(results, expected_groups):
    assert results.groups == expected_groups


def test_rename_tax_ids(samples):
    taxa_df = pd.DataFrame(
        {"tax1": [10, 20], "tax2": [30, 40], "unknown": [50, 60]}, index=["s1", "s2"]
    )

    with mock.patch.object(
        type(samples),
        "taxonomy",
        new_callable=mock.PropertyMock,
        return_value={"name": {"tax1": "Bacteroides", "tax2": "Prevotella"}},
    ):
        result = samples._rename_tax_ids(taxa_df)

    assert list(result.columns) == ["Bacteroides (tax1)", "Prevotella (tax2)", "unknown"]
    assert result.values.tolist() == taxa_df.values.tolist()


def test_drop_samples_with_zero_abundance(samples):
    metadata_df = pd.DataFrame({"group": ["a", "b", "c"]}, index=["s1", "s2", "s3"])
    taxa_df = pd.DataFrame({"t1": [1, 0, 3], "t2": [4, 0, 6]}, index=["s1", "s2", "s3"])

    with pytest.warns(StatsWarning, match="1 sample was excluded.*zero abundance"):
        result = samples._drop_samples_with_zero_abundance(metadata_df, taxa_df)

    assert list(result.index) == ["s1", "s3"]


def test_drop_samples_with_zero_abundance_none_dropped(samples):
    metadata_df = pd.DataFrame({"group": ["a", "b"]}, index=["s1", "s2"])
    taxa_df = pd.DataFrame({"t1": [1, 2], "t2": [3, 4]}, index=["s1", "s2"])

    result = samples._drop_samples_with_zero_abundance(metadata_df, taxa_df)

    assert list(result.index) == ["s1", "s2"]


def test_assert_min_num_groups_for_global_test(samples):
    df = pd.DataFrame({"group": ["a", "b", "a", "b"]})

    with pytest.raises(StatsException, match="include_global_test.*at least 3 groups"):
        samples._assert_min_num_groups_for_global_test(df, "group")


def _assert_ancombc_covariates(
    main_results, group_by_variable, reference_group, non_reference_groups
):
    covariates = set(main_results.index.get_level_values("Covariate").unique())
    expected_covariates = {"Intercept"} | {
        f"{group_by_variable}: {g} vs {reference_group} (reference)" for g in non_reference_groups
    }
    assert covariates == expected_covariates


def test_ancombc(ocx, api_data, samples):
    samples.extend([ocx.Samples.get("cc18208d98ad48b3"), ocx.Samples.get("5445740666134eee")])

    for sample, group in zip(samples, ["a", "b", "a", "b", "a"]):
        sample.metadata.custom["group"] = group

    results = samples._ancombc(group_by="group", metric="readcount_w_children", rank="genus")

    assert isinstance(results, AncombcResults)
    assert results.reference_group == "a"
    assert results.alpha == 0.05
    assert results.adjustment_method == AdjustmentMethod.BenjaminiHochberg
    assert results.sample_size == 5
    assert results.group_by_variable == "group"
    assert results.group_sizes == {"a": 3, "b": 2}
    assert results.global_results is None
    assert set(results.main_results.columns) == ANCOMBC_MAIN_RESULTS_COLUMNS
    assert results.main_results.index.names == ANCOMBC_MAIN_RESULTS_INDEX_NAMES
    assert len(results.main_results) > 0
    _assert_ancombc_covariates(results.main_results, "group", "a", ["b"])


def test_ancombc_composite_group_by(ocx, api_data, samples):
    samples.extend([ocx.Samples.get("cc18208d98ad48b3"), ocx.Samples.get("5445740666134eee")])

    for sample, c1, c2 in zip(samples, ["x", "x", "y", "y", "x"], ["1", "1", "2", "2", "1"]):
        sample.metadata.custom["col1"] = c1
        sample.metadata.custom["col2"] = c2

    results = samples._ancombc(
        group_by=("col1", "col2"), metric="readcount_w_children", rank="genus"
    )

    assert results.group_by_variable == "col1_col2"
    assert results.group_sizes == {"x_1": 3, "y_2": 2}
    assert results.reference_group == "x_1"
    assert set(results.main_results.columns) == ANCOMBC_MAIN_RESULTS_COLUMNS
    assert results.main_results.index.names == ANCOMBC_MAIN_RESULTS_INDEX_NAMES
    _assert_ancombc_covariates(results.main_results, "col1_col2", "x_1", ["y_2"])


def test_ancombc_explicit_reference_group(ocx, api_data, samples):
    samples.extend([ocx.Samples.get("cc18208d98ad48b3"), ocx.Samples.get("5445740666134eee")])

    for sample, group in zip(samples, ["a", "b", "a", "b", "a"]):
        sample.metadata.custom["group"] = group

    results = samples._ancombc(
        group_by="group",
        reference_group="b",
        metric="readcount_w_children",
        rank="genus",
    )

    assert results.reference_group == "b"
    _assert_ancombc_covariates(results.main_results, "group", "b", ["a"])


def test_ancombc_invalid_reference_group(ocx, api_data, samples):
    samples.extend([ocx.Samples.get("cc18208d98ad48b3"), ocx.Samples.get("5445740666134eee")])

    for sample, group in zip(samples, ["a", "b", "a", "b", "a"]):
        sample.metadata.custom["group"] = group

    with pytest.raises(StatsException, match="reference_group.*group name"):
        samples._ancombc(
            group_by="group",
            reference_group="nonexistent",
            metric="readcount_w_children",
            rank="genus",
        )


def test_ancombc_with_global_test(ocx, api_data, samples):
    samples.extend(
        [
            ocx.Samples.get("cc18208d98ad48b3"),
            ocx.Samples.get("5445740666134eee"),
            ocx.Samples.get("0ecac25ec0004fe4"),
        ]
    )

    for sample, group in zip(samples, ["a", "b", "c", "a", "b", "c"]):
        sample.metadata.custom["group"] = group

    with pytest.warns(UserWarning):
        results = samples._ancombc(
            group_by="group",
            include_global_test=True,
            metric="readcount_w_children",
            rank="genus",
            require_classification_version_match=False,
        )

    assert results.reference_group == "a"
    assert results.group_sizes == {"a": 2, "b": 2, "c": 2}
    assert set(results.main_results.columns) == ANCOMBC_MAIN_RESULTS_COLUMNS
    assert results.main_results.index.names == ANCOMBC_MAIN_RESULTS_INDEX_NAMES
    _assert_ancombc_covariates(results.main_results, "group", "a", ["b", "c"])
    assert results.global_results is not None
    assert results.global_results.index.name == ANCOMBC_GLOBAL_RESULTS_INDEX_NAME
    assert set(results.global_results.columns) == ANCOMBC_GLOBAL_RESULTS_COLUMNS


def test_ancombc_global_test_too_few_groups(ocx, api_data, samples):
    samples.extend([ocx.Samples.get("cc18208d98ad48b3"), ocx.Samples.get("5445740666134eee")])

    for sample, group in zip(samples, ["a", "b", "a", "b", "a"]):
        sample.metadata.custom["group"] = group

    with pytest.raises(StatsException, match="include_global_test.*at least 3 groups"):
        samples._ancombc(
            group_by="group",
            include_global_test=True,
            metric="readcount_w_children",
            rank="genus",
        )


def test_ancombc_three_groups(ocx, api_data, samples):
    samples.extend(
        [
            ocx.Samples.get("cc18208d98ad48b3"),
            ocx.Samples.get("5445740666134eee"),
            ocx.Samples.get("0ecac25ec0004fe4"),
        ]
    )

    for sample, group in zip(samples, ["a", "b", "c", "a", "b", "c"]):
        sample.metadata.custom["group"] = group

    with pytest.warns(UserWarning):
        results = samples._ancombc(
            group_by="group",
            metric="readcount_w_children",
            rank="genus",
            require_classification_version_match=False,
        )

    assert results.reference_group == "a"
    assert results.sample_size == 6
    assert results.group_sizes == {"a": 2, "b": 2, "c": 2}
    assert set(results.main_results.columns) == ANCOMBC_MAIN_RESULTS_COLUMNS
    assert results.main_results.index.names == ANCOMBC_MAIN_RESULTS_INDEX_NAMES
    _assert_ancombc_covariates(results.main_results, "group", "a", ["b", "c"])
