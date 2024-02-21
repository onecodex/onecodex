from __future__ import annotations
from dataclasses import dataclass
import warnings
from typing import TYPE_CHECKING, Optional

from onecodex.exceptions import OneCodexException, StatsException, StatsWarning
from onecodex.lib.enums import Rank, AlphaDiversityMetric, AlphaDiversityStatsTest

if TYPE_CHECKING:
    import pandas as pd


@dataclass(frozen=True)
class AlphaDiversityStatsResults:
    test: AlphaDiversityStatsTest
    statistic: float
    pvalue: float
    group_by_variable: str
    groups: set[str]
    paired_by_variable: Optional[str] = None
    posthoc_df: Optional[pd.DataFrame] = None


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
            Metadata variable to pair samples in each group. May only be used with `test="wilcoxon"`.
            If `paired_by` is a tuple or list, field values are joined with an underscore character
            ("_").
        test : {'auto', 'wilcoxon', 'mannwhitneyu', 'kruskal'}, optional
            Stats test to perform. If `'auto'`, `'mannwhitneyu'` will be chosen if there are two
            groups of unpaired data. `'wilcoxon'` will be chosen if there are two groups and
            `paired_by` is specified. `'kruskal'` will be chosen if there are more than 2 groups.
        metric : {'shannon', 'simpson', 'observed_taxa'}, optional
            The alpha diversity metric to calculate.
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.
        alpha : float, optional
            Threshold to determine statistical significance when `test="kruskal"` (e.g. p < `alpha`).
            Must be between 0 and 1 (exclusive). If the Kruskal-Wallis p-value is significant and
            there are more than two groups, a posthoc Dunn test is performed.

        Returns
        -------
        AlphaDiversityStatsResults
            A dataclass with these attributes:
            - `test`: the stats test that was performed
            - `statistic`: the computed test statistic (e.g. U statistic if `test="mannwhitneyu"`)
            - `pvalue`: the computed p-value
            - `group_by_variable`: the name of the variable used to group samples by
            - `groups`: the names of the groups defined by `group_by_variable`
            - `paired_by_variable`: the name of the variable used to pair samples by (if the data
              were paired)
            - `posthoc_df`: `pd.DataFrame` containing Dunn test p-values (if a Dunn test was
              performed). p-values are adjusted for false discovery rate using Benjamini-Hochberg.
              The index and columns are sorted group names.

        See Also
        --------
        scipy.stats.wilcoxon
        scipy.stats.mannwhitneyu
        scipy.stats.kruskal
        scikit_posthocs.posthoc_dunn

        """
        if not (0.0 < alpha < 1.0):
            raise StatsException("`alpha` must be between 0 and 1 (exclusive).")

        group_by = self._tuplize(group_by)
        metadata_fields = [group_by]
        if paired_by is not None:
            paired_by = self._tuplize(paired_by)
            metadata_fields.append(paired_by)

        # Munge the metadata first in case there's any errors
        df, magic_fields = self._metadata_fetch(
            metadata_fields, coerce_missing_composite_fields=False
        )
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
        prev_count = len(df)
        df.dropna(subset=[group_by_column_name], inplace=True)
        num_dropped = prev_count - len(df)
        if num_dropped > 0:
            warnings.warn(
                f"{num_dropped} sample{'s were' if num_dropped > 1 else ' was'} excluded from the "
                f"test as {'they were' if num_dropped > 1 else 'it was'} missing `group_by` "
                f"metadata.",
                StatsWarning,
            )

        if paired_by_column_name:
            # Drop samples with missing paired_by data
            prev_count = len(df)
            df.dropna(subset=[paired_by_column_name], inplace=True)
            num_dropped = prev_count - len(df)
            if num_dropped > 0:
                warnings.warn(
                    f"{num_dropped} sample{'s were' if num_dropped > 1 else ' was'} excluded from "
                    f"the test as {'they were' if num_dropped > 1 else 'it was'} missing "
                    f"`paired_by` metadata.",
                    StatsWarning,
                )

        # Drop groups of size < 2
        prev_count = len(df)
        df = df.groupby(group_by_column_name).filter(lambda group: len(group) > 1)
        num_dropped = prev_count - len(df)
        if num_dropped > 0:
            warnings.warn(
                f"{num_dropped} sample{'s were' if num_dropped > 1 else ' was'} excluded from the "
                f"test as {'they' if num_dropped > 1 else 'it'} belonged to group(s) of size < 2.",
                StatsWarning,
            )

        num_groups = df[group_by_column_name].nunique()
        if num_groups < 2:
            raise StatsException(
                f"`group_by` must have at least two groups to test between after filtering (found {num_groups})."
            )

        if test == AlphaDiversityStatsTest.Auto:
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

        self._assert_two_groups(df, group_by_column_name)

        (group1_name, group1_df), (group2_name, group2_df) = df.groupby(group_by_column_name)
        if len(group1_df) != len(group2_df):
            raise StatsException(
                f"`group_by` must have two groups of the same size after filtering (found size {len(group1_df)} and size {len(group2_df)})."
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
            group_by_variable=group_by_column_name,
            groups={group1_name, group2_name},
            paired_by_variable=paired_by_column_name,
        )

    def _mannwhitneyu(
        self, df: pd.DataFrame, metric: str, group_by_column_name: str
    ) -> AlphaDiversityStatsResults:
        from scipy.stats import mannwhitneyu

        self._assert_two_groups(df, group_by_column_name)

        group_names, group_alpha_values = self._get_group_names_and_alpha_values(
            df, metric, group_by_column_name
        )
        result = mannwhitneyu(*group_alpha_values)

        return AlphaDiversityStatsResults(
            test=AlphaDiversityStatsTest.Mannwhitneyu,
            statistic=result.statistic,
            pvalue=result.pvalue,
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

        posthoc_df = None
        if result.pvalue < alpha and len(group_names) > 2:
            posthoc_df = posthoc_dunn(
                df, val_col=metric, group_col=group_by_column_name, p_adjust="fdr_bh"
            )

        return AlphaDiversityStatsResults(
            test=AlphaDiversityStatsTest.Kruskal,
            statistic=result.statistic,
            pvalue=result.pvalue,
            group_by_variable=group_by_column_name,
            groups=group_names,
            posthoc_df=posthoc_df,
        )

    def _tuplize(self, value) -> tuple:
        if not isinstance(value, tuple):
            if isinstance(value, list):
                # Convert to tuple so that it's hashable in self._metadata_fetch()
                value = tuple(value)
            else:
                value = (value,)
        return value

    def _assert_two_groups(self, df: pd.DataFrame, group_by_column_name: str):
        num_groups = df[group_by_column_name].nunique()
        if num_groups != 2:
            raise StatsException(
                f"`group_by` must have exactly two groups to test between after filtering (found {num_groups})."
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
