from __future__ import annotations
import itertools
import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

from onecodex.exceptions import OneCodexException, StatsException, StatsWarning
from onecodex.lib.enums import (
    Rank,
    AlphaDiversityMetric,
    AlphaDiversityStatsTest,
    BetaDiversityMetric,
    BetaDiversityStatsTest,
)

if TYPE_CHECKING:
    import pandas as pd
    import skbio


# TODO when min Python version is >=3.10, create a base class with
# @dataclass(frozen=True, kw_only=True) and have AlphaDiversityStatsResults and
# BetaDiversityStatsResults inherit from it


@dataclass(frozen=True)
class AlphaDiversityStatsResults:
    test: AlphaDiversityStatsTest
    statistic: float
    pvalue: float
    sample_size: int
    group_by_variable: str
    groups: set[str]
    paired_by_variable: Optional[str] = None
    posthoc: Optional[PosthocResults] = None

    @property
    def posthoc_df(self) -> Optional[pd.DataFrame]:
        warnings.warn(
            "`posthoc_df` is deprecated and will be removed in a future release. Please use "
            "`posthoc.adjusted_pvalues` instead.",
            DeprecationWarning,
        )
        return self.posthoc.adjusted_pvalues if self.posthoc else None


@dataclass(frozen=True)
class BetaDiversityStatsResults:
    test: BetaDiversityStatsTest
    statistic: float
    pvalue: float
    num_permutations: int
    sample_size: int
    group_by_variable: str
    groups: set[str]
    posthoc: Optional[PosthocResults] = None


@dataclass(frozen=True)
class PosthocResults:
    adjusted_pvalues: pd.DataFrame
    pvalues: Optional[pd.DataFrame] = None
    statistics: Optional[pd.DataFrame] = None


class StatsMixin:
    def alpha_diversity_stats(
        self,
        *,
        group_by: str | tuple[str, ...] | list[str],
        paired_by: Optional[str | tuple[str, ...] | list[str]] = None,
        test: AlphaDiversityStatsTest = AlphaDiversityStatsTest.Auto,
        metric: AlphaDiversityMetric = AlphaDiversityMetric.Shannon,
        rank: Rank = Rank.Auto,
        alpha: float = 0.05,
    ) -> AlphaDiversityStatsResults:
        """Perform a test for significant differences between groups of alpha diversity values.

        The following tests are supported:

        - Wilcoxon (2 groups, paired data)
        - Mann-Whitney U (2 groups, unpaired data)
        - Kruskal-Wallis with optional posthoc Dunn test (>=2 groups, unpaired data)

        Parameters
        ----------
        group_by : str or tuple of str or list of str
            Metadata variable to group samples by. At least two groups are required. If `group_by`
            is a tuple or list, field values are joined with an underscore character ("_").
        paired_by : str or tuple of str or list of str, optional
            Metadata variable to pair samples in each group. May only be used with
            `test="wilcoxon"`. If `paired_by` is a tuple or list, field values are joined with an
            underscore character ("_").
        test : {'auto', 'wilcoxon', 'mannwhitneyu', 'kruskal'}, optional
            Stats test to perform. If `'auto'`, `'mannwhitneyu'` will be chosen if there are two
            groups of unpaired data. `'wilcoxon'` will be chosen if there are two groups and
            `paired_by` is specified. `'kruskal'` will be chosen if there are more than 2 groups.
        metric : {'shannon', 'simpson', 'observed_taxa'}, optional
            The alpha diversity metric to calculate. Note that Shannon diversity is calculated using
            log base ``e`` (natural log).
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.
        alpha : float, optional
            Threshold to determine statistical significance when `test="kruskal"`
            (e.g. p < `alpha`). Must be between 0 and 1 (exclusive). If the Kruskal-Wallis p-value
            is significant and there are more than two groups, a posthoc Dunn test is performed.

        Returns
        -------
        AlphaDiversityStatsResults
            A dataclass with these attributes:
            - `test`: stats test that was performed
            - `statistic`: computed test statistic (e.g. U statistic if `test="mannwhitneyu"`)
            - `pvalue`: computed p-value
            - `sample_size`: number of samples used in the test after filtering
            - `group_by_variable`: name of the variable used to group samples by
            - `groups`: names of the groups defined by `group_by_variable`
            - `paired_by_variable`: name of the variable used to pair samples by (if the data were
              paired)
            - `posthoc`: A dataclass with these attributes (if posthoc results were computed):
              - `adjusted_pvalues`: `pd.DataFrame` containing Dunn test *adjusted* p-values.
                p-values are adjusted for false discovery rate using Benjamini-Hochberg. The index
                and columns are sorted group names.

        See Also
        --------
        scipy.stats.wilcoxon
        scipy.stats.mannwhitneyu
        scipy.stats.kruskal
        scikit_posthocs.posthoc_dunn

        """
        self._assert_valid_alpha(alpha)

        group_by = self._tuplize(group_by)
        metadata_fields = [group_by]
        if paired_by is not None:
            paired_by = self._tuplize(paired_by)
            metadata_fields.append(paired_by)

        # Munge the metadata first in case there's any errors
        metadata_results = self._metadata_fetch(
            metadata_fields, coerce_missing_composite_fields=False
        )
        df = metadata_results.df
        magic_fields = metadata_results.renamed_fields
        group_by_column_name = magic_fields[group_by]
        paired_by_column_name = magic_fields.get(paired_by, None)

        # Compute alpha diversity
        df.loc[:, metric] = self.alpha_diversity(metric=metric, rank=rank)

        # Drop NaN alpha diversity values
        prev_count = len(df)
        df.dropna(subset=[metric], inplace=True)
        num_dropped = prev_count - len(df)
        if num_dropped > 0:
            warnings.warn(
                f"{num_dropped} sample{'s were' if num_dropped > 1 else ' was'} excluded from the "
                f"test as no {metric} value was calculated for "
                f"{'them' if num_dropped > 1 else 'it'}.",
                StatsWarning,
            )

        # Drop samples with missing group_by data
        df = self._drop_missing_data(df, group_by_column_name, "group_by")

        if paired_by_column_name:
            # Drop samples with missing paired_by data
            df = self._drop_missing_data(df, paired_by_column_name, "paired_by")

        # Drop samples missing abundance data
        df = self._drop_classifications_without_abundances(df)

        # Drop groups of size < 2
        # NOTE: it's important to do this filtering *after* the other filters in case those filters
        # change the group sizes.
        df = self._drop_group_sizes_smaller_than(df, group_by_column_name, 2)

        self._assert_min_num_groups(df, group_by_column_name, 2)

        if test == AlphaDiversityStatsTest.Auto:
            num_groups = df[group_by_column_name].nunique()
            if num_groups == 2:
                if paired_by_column_name:
                    test = AlphaDiversityStatsTest.Wilcoxon
                else:
                    test = AlphaDiversityStatsTest.Mannwhitneyu
            else:
                test = AlphaDiversityStatsTest.Kruskal

        if paired_by_column_name and test != AlphaDiversityStatsTest.Wilcoxon:
            raise StatsException('`paired_by` may only be specified with `test="wilcoxon"`.')

        if test == AlphaDiversityStatsTest.Wilcoxon:
            return self._wilcoxon(df, metric, group_by_column_name, paired_by_column_name)
        elif test == AlphaDiversityStatsTest.Mannwhitneyu:
            return self._mannwhitneyu(df, metric, group_by_column_name)
        elif test == AlphaDiversityStatsTest.Kruskal:
            return self._kruskal(df, metric, group_by_column_name, alpha)
        else:
            raise OneCodexException(f"Stats test {test} is not supported.")

    def _assert_valid_alpha(self, alpha: float):
        if not (0.0 < alpha < 1.0):
            raise StatsException("`alpha` must be between 0 and 1 (exclusive).")

    def _tuplize(self, value) -> tuple:
        if not isinstance(value, tuple):
            if isinstance(value, list):
                # Convert to tuple so that it's hashable in self._metadata_fetch()
                value = tuple(value)
            else:
                value = (value,)
        return value

    def _drop_missing_data(
        self, df: pd.DataFrame, column_name: str, param_name: str
    ) -> pd.DataFrame:
        prev_count = len(df)
        df = df.dropna(subset=[column_name])
        num_dropped = prev_count - len(df)
        if num_dropped > 0:
            warnings.warn(
                f"{num_dropped} sample{'s were' if num_dropped > 1 else ' was'} excluded from the "
                f"test as {'they were' if num_dropped > 1 else 'it was'} missing `{param_name}` "
                f"metadata.",
                StatsWarning,
            )
        return df

    def _drop_classifications_without_abundances(self, df: pd.DataFrame) -> pd.DataFrame:
        prev_count = len(df)
        df = df.drop(
            [id_ for id_ in self._classification_ids_without_abundances if id_ in df.index]
        )
        num_dropped = prev_count - len(df)
        if num_dropped > 0:
            warnings.warn(
                f"{num_dropped} sample{'s were' if num_dropped > 1 else ' was'} excluded from the "
                f"test as {'they were' if num_dropped > 1 else 'it was'} missing abundance "
                f"calculations.",
                StatsWarning,
            )
        return df

    def _drop_group_sizes_smaller_than(
        self, df: pd.DataFrame, group_by_column_name: str, min_group_size: int
    ) -> pd.DataFrame:
        prev_count = len(df)
        df = df.groupby(group_by_column_name).filter(lambda group: len(group) >= min_group_size)
        num_dropped = prev_count - len(df)
        if num_dropped > 0:
            warnings.warn(
                f"{num_dropped} sample{'s were' if num_dropped > 1 else ' was'} excluded from the "
                f"test as {'they' if num_dropped > 1 else 'it'} belonged to group(s) of size < "
                f"{min_group_size}.",
                StatsWarning,
            )
        return df

    def _assert_min_num_groups(
        self, df: pd.DataFrame, group_by_column_name: str, min_num_groups: int
    ):
        num_groups = df[group_by_column_name].nunique()
        if num_groups < min_num_groups:
            raise StatsException(
                f"`group_by` must have at least {min_num_groups} groups to test between after "
                f"filtering (found {num_groups})."
            )

    def _assert_exact_num_groups(
        self, df: pd.DataFrame, group_by_column_name: str, exact_num_groups: int
    ):
        num_groups = df[group_by_column_name].nunique()
        if num_groups != exact_num_groups:
            raise StatsException(
                f"`group_by` must have exactly {exact_num_groups} groups to test between after "
                f"filtering (found {num_groups})."
            )

    def _get_group_names_and_alpha_values(
        self, df: pd.DataFrame, metric: str, group_by_column_name: str
    ) -> tuple[set[str], list[pd.Series]]:
        group_names = set()
        group_alpha_values = []
        for group_name, group_df in df.groupby(group_by_column_name):
            group_names.add(group_name)
            group_alpha_values.append(group_df[metric])

        return group_names, group_alpha_values

    def _wilcoxon(
        self,
        df: pd.DataFrame,
        metric: str,
        group_by_column_name: str,
        paired_by_column_name: Optional[str],
    ) -> AlphaDiversityStatsResults:
        from scipy.stats import wilcoxon

        if not paired_by_column_name:
            raise StatsException('`paired_by` must be specified with `test="wilcoxon"`.')

        self._assert_exact_num_groups(df, group_by_column_name, 2)

        (group1_name, group1_df), (group2_name, group2_df) = df.groupby(group_by_column_name)
        if len(group1_df) != len(group2_df):
            raise StatsException(
                f"`group_by` must have two groups of the same size after filtering "
                f"(found size {len(group1_df)} and size {len(group2_df)})."
            )

        if group1_df[paired_by_column_name].nunique() != len(group1_df) or group2_df[
            paired_by_column_name
        ].nunique() != len(group2_df):
            raise StatsException("`paired_by` must have unique values within each group.")

        if set(group1_df[paired_by_column_name]) != set(group2_df[paired_by_column_name]):
            raise StatsException(
                "Every sample in one group must have a matching paired sample in the other group."
            )

        # Reorder one of the group's samples based on the ordering of the `paired_by` column in the
        # other group. The actual ordering does not matter.
        group2_df.set_index(paired_by_column_name, inplace=True)
        group2_df = group2_df.reindex(index=group1_df[paired_by_column_name])

        result = wilcoxon(group1_df[metric], group2_df[metric])

        return AlphaDiversityStatsResults(
            test=AlphaDiversityStatsTest.Wilcoxon,
            statistic=result.statistic,
            pvalue=result.pvalue,
            sample_size=len(df),
            group_by_variable=group_by_column_name,
            groups={group1_name, group2_name},
            paired_by_variable=paired_by_column_name,
        )

    def _mannwhitneyu(
        self, df: pd.DataFrame, metric: str, group_by_column_name: str
    ) -> AlphaDiversityStatsResults:
        from scipy.stats import mannwhitneyu

        self._assert_exact_num_groups(df, group_by_column_name, 2)

        group_names, group_alpha_values = self._get_group_names_and_alpha_values(
            df, metric, group_by_column_name
        )
        result = mannwhitneyu(*group_alpha_values)

        return AlphaDiversityStatsResults(
            test=AlphaDiversityStatsTest.Mannwhitneyu,
            statistic=result.statistic,
            pvalue=result.pvalue,
            sample_size=len(df),
            group_by_variable=group_by_column_name,
            groups=group_names,
        )

    def _kruskal(
        self, df: pd.DataFrame, metric: str, group_by_column_name: str, alpha: float
    ) -> AlphaDiversityStatsResults:
        from scikit_posthocs import posthoc_dunn
        from scipy.stats import kruskal

        group_names, group_alpha_values = self._get_group_names_and_alpha_values(
            df, metric, group_by_column_name
        )
        result = kruskal(*group_alpha_values)

        posthoc = None
        if result.pvalue < alpha and len(group_names) > 2:
            posthoc = PosthocResults(
                adjusted_pvalues=posthoc_dunn(
                    df, val_col=metric, group_col=group_by_column_name, p_adjust="fdr_bh"
                )
            )

        return AlphaDiversityStatsResults(
            test=AlphaDiversityStatsTest.Kruskal,
            statistic=result.statistic,
            pvalue=result.pvalue,
            sample_size=len(df),
            group_by_variable=group_by_column_name,
            groups=group_names,
            posthoc=posthoc,
        )

    def beta_diversity_stats(
        self,
        *,
        group_by: str | tuple[str, ...] | list[str],
        metric: BetaDiversityMetric = BetaDiversityMetric.BrayCurtis,
        rank: Rank = Rank.Auto,
        alpha: float = 0.05,
        num_permutations: int = 999,
    ) -> BetaDiversityStatsResults:
        """Test for significant differences between groups of samples based on their distances.

        Beta diversity distances between samples are computed and a PERMANOVA test is performed to
        assess whether there are significant differences between groups of samples. Posthoc pairwise
        PERMANOVA tests are performed if the global test is found to be statistically significant
        and there are more than two groups.

        Parameters
        ----------
        group_by : str or tuple of str or list of str
            Metadata variable to group samples by. At least two groups are required. If `group_by`
            is a tuple or list, field values are joined with an underscore character ("_").
        metric : {'braycurtis', 'jaccard', 'cityblock', 'manhattan', 'aitchison', 'unweighted_unifrac', 'weighted_unifrac'}, optional
            The beta diversity distance metric to calculate. Note that 'cityblock' and 'manhattan'
            are equivalent metrics.
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.
        alpha : float, optional
            Threshold to determine statistical significance (e.g. p < `alpha`). Must be between 0
            and 1 (exclusive). If the p-value is significant and there are more than two groups,
            posthoc pairwise PERMANOVA tests are performed.
        num_permutations : int, optional
            Number of permutations to use when computing the p-value.

        Returns
        -------
        BetaDiversityStatsResults
            A dataclass with these attributes:
            - `test`: stats test that was performed
            - `statistic`: PERMANOVA pseudo-F test statistic
            - `pvalue`: p-value based on `num_permutations`
            - `num_permutations`: number of permutations used to compute `pvalue`
            - `sample_size`: number of samples used in the test after filtering
            - `group_by_variable`: name of the variable used to group samples by
            - `groups`: names of the groups defined by `group_by_variable`
            - `posthoc`: A dataclass with these attributes (if posthoc results were computed):
              - `statistics`: `pd.DataFrame` containing pairwise PERMANOVA pseudo-F statistics. The
                index and columns are sorted group names.
              - `pvalues`: `pd.DataFrame` containing pairwise PERMANOVA *unadjusted* p-values. The
                index and columns are sorted group names.
              - `adjusted_pvalues`: `pd.DataFrame` containing pairwise PERMANOVA *adjusted*
                p-values. p-values are adjusted for false discovery rate using Benjamini-Hochberg.
                The index columns are sorted group names.

        See Also
        --------
        skbio.stats.distance.permanova
        scipy.stats.false_discovery_control

        """
        self._assert_valid_alpha(alpha)

        # Munge the metadata first in case there's any errors
        group_by = self._tuplize(group_by)
        metadata_results = self._metadata_fetch([group_by], coerce_missing_composite_fields=False)
        df = metadata_results.df
        magic_fields = metadata_results.renamed_fields
        group_by_column_name = magic_fields[group_by]

        # Drop samples with missing group_by data
        df = self._drop_missing_data(df, group_by_column_name, "group_by")

        # Drop samples missing abundance data
        df = self._drop_classifications_without_abundances(df)

        # Drop groups of size < 2
        # NOTE: it's important to do this filtering *after* the other filters in case those filters
        # change the group sizes.
        df = self._drop_group_sizes_smaller_than(df, group_by_column_name, 2)

        self._assert_min_num_groups(df, group_by_column_name, 2)

        # Compute beta diversity distance matrix and filter it to match the remaining IDs in the
        # metadata dataframe
        dm = self.beta_diversity(metric=metric, rank=rank).filter(df.index, strict=True)

        return self._permanova(dm, df, group_by_column_name, alpha, num_permutations)

    def _permanova(
        self,
        dm: skbio.DistanceMatrix,
        df: pd.DataFrame,
        group_by_column_name: str,
        alpha: float,
        num_permutations: int,
    ) -> BetaDiversityStatsResults:
        import numpy as np
        import pandas as pd
        from scipy.stats import false_discovery_control
        from skbio.stats.distance import permanova

        group_names = {name for name, _ in df.groupby(group_by_column_name)}
        result = permanova(dm, df, column=group_by_column_name, permutations=num_permutations)

        posthoc = None
        if result["p-value"] < alpha and result["number of groups"] > 2:
            # Compute PERMANOVA for each pair of groups. Assuming we have groups "g1", "g2", "g3",
            # we'll run PERMANOVA 3 times: g1 vs g2, g1 vs g3, g2 vs g3. We'll store the results in
            # a 2D matrix ordered by sorted group names. Since results are symmetric and we don't
            # need to compute the diagonal (e.g. g1 vs g1), we can just compute the upper triangle
            # in row-major order:
            #
            #   g1  g2  g3
            # g1 x   1   2
            # g2 x   x   3
            # g3 x   x   x
            group_names = sorted(group_names)
            upper_triangle_statistics = []
            upper_triangle_pvalues = []
            for group1, group2 in itertools.combinations(group_names, 2):
                pairwise_dm = dm.filter(
                    df[
                        (df[group_by_column_name] == group1) | (df[group_by_column_name] == group2)
                    ].index,
                    strict=True,
                )
                # The dataframe can be a superset of the distance matrix and the IDs don't need
                # to be in the same order.
                pairwise_result = permanova(
                    pairwise_dm, df, column=group_by_column_name, permutations=num_permutations
                )
                assert pairwise_result["number of groups"] == 2

                upper_triangle_statistics.append(pairwise_result["test statistic"])
                upper_triangle_pvalues.append(pairwise_result["p-value"])

            num_groups = len(group_names)
            shape = (num_groups, num_groups)
            upper_triangle_indices = np.triu_indices(num_groups, k=1)

            statistics = np.zeros(shape, dtype=float)
            statistics[upper_triangle_indices] = upper_triangle_statistics
            statistics += statistics.T
            np.fill_diagonal(statistics, np.nan)

            pvalues = np.zeros(shape, dtype=float)
            pvalues[upper_triangle_indices] = upper_triangle_pvalues
            pvalues += pvalues.T
            np.fill_diagonal(pvalues, 1.0)

            adjusted_pvalues = np.zeros(shape, dtype=float)
            adjusted_pvalues[upper_triangle_indices] = false_discovery_control(
                upper_triangle_pvalues, method="bh"
            )
            adjusted_pvalues += adjusted_pvalues.T
            np.fill_diagonal(adjusted_pvalues, 1.0)

            posthoc = PosthocResults(
                statistics=pd.DataFrame(statistics, index=group_names, columns=group_names),
                pvalues=pd.DataFrame(pvalues, index=group_names, columns=group_names),
                adjusted_pvalues=pd.DataFrame(
                    adjusted_pvalues, index=group_names, columns=group_names
                ),
            )

        return BetaDiversityStatsResults(
            test=BetaDiversityStatsTest.Permanova,
            statistic=result["test statistic"],
            pvalue=result["p-value"],
            num_permutations=result["number of permutations"],
            sample_size=result["sample size"],
            group_by_variable=group_by_column_name,
            groups=group_names,
            posthoc=posthoc,
        )
