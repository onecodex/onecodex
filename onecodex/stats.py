from __future__ import annotations
import itertools
import re
import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

from onecodex.distance import DistanceMixin
from onecodex.exceptions import OneCodexException, StatsException, StatsWarning
from onecodex.lib.enums import (
    Metric,
    Rank,
    AlphaDiversityMetric,
    AlphaDiversityStatsTest,
    BetaDiversityMetric,
    BetaDiversityStatsTest,
    AdjustmentMethod,
    PosthocStatsTest,
)
from onecodex.models.base_sample_collection import BaseSampleCollection

if TYPE_CHECKING:
    import pandas as pd
    import skbio


@dataclass(frozen=True, kw_only=True)
class StatsResults:
    alpha: float
    sample_size: int
    group_by_variable: str
    group_sizes: dict[str, int]

    @property
    def groups(self) -> set[str]:
        return set(self.group_sizes)


@dataclass(frozen=True, kw_only=True)
class DiversityStatsResults(StatsResults):
    statistic: float
    pvalue: float
    posthoc: Optional[PosthocResults] = None


@dataclass(frozen=True, kw_only=True)
class AlphaDiversityStatsResults(DiversityStatsResults):
    """A dataclass for storing the results of an alpha diversity stats test.

    - `test`: stats test that was performed
    - `statistic`: computed test statistic (e.g. U statistic if `test="mannwhitneyu"`)
    - `pvalue`: computed p-value
    - `alpha`: p-value threshold used to determine whether to run a posthoc test when `test="kruskal"`
    - `sample_size`: number of samples used in the test after filtering
    - `group_by_variable`: name of the variable used to group samples by
    - `group_sizes`: dict mapping group name to sample size in each group
    - `paired_by_variable`: name of the variable used to pair samples by (if the data were
    paired)
    - `posthoc`: :class:`~onecodex.stats.PosthocResults`
    """

    test: AlphaDiversityStatsTest
    paired_by_variable: Optional[str] = None


@dataclass(frozen=True, kw_only=True)
class BetaDiversityStatsResults(DiversityStatsResults):
    """A dataclass for storing the results of a beta diversity test.

    - `test`: stats test that was performed
    - `statistic`: PERMANOVA pseudo-F test statistic
    - `pvalue`: p-value based on `num_permutations`
    - `alpha`: p-value threshold used to determine whether to run a posthoc test
    - `num_permutations`: number of permutations used to compute `pvalue`
    - `sample_size`: number of samples used in the test after filtering
    - `group_by_variable`: name of the variable used to group samples by
    - `group_sizes`: dict mapping group name to sample size in each group
    - `posthoc`: :class:`~onecodex.stats.PosthocResults`
    """

    test: BetaDiversityStatsTest
    num_permutations: int


@dataclass(frozen=True, kw_only=True)
class PosthocResults:
    """A dataclass for storing results from a posthoc statistical test.

    - `test`: posthoc stats test that was performed
    - `adjustment_method`: method used to adjust p-values to control the false discovery rate
    - `statistics`: `pd.DataFrame` containing pairwise test statistics. The index and columns are
    sorted group names.
    - `pvalues`: `pd.DataFrame` containing *unadjusted* p-values. The index and columns are sorted
    group names.
    - `adjusted_pvalues`: `pd.DataFrame` containing *adjusted* p-values. p-values are adjusted for
    false discovery rate using `adjustment_method`. The index and columns are sorted group names.
    """

    test: PosthocStatsTest
    adjustment_method: AdjustmentMethod
    adjusted_pvalues: pd.DataFrame
    pvalues: Optional[pd.DataFrame] = None
    statistics: Optional[pd.DataFrame] = None


@dataclass(frozen=True, kw_only=True)
class AncombcResults(StatsResults):
    """A dataclass for storing the results of an ANCOM-BC differential abundance test.

    - `main_results`: `pd.DataFrame` containing ANCOM-BC test results for each taxon, including
    log2-fold change (with standard error), W-statistic, p-value, adjusted p-value, and whether or
    not the result is statistically significant. See `skbio.stats.composition.ancombc` for details.
    - `global_results`: `pd.DataFrame` containing global ANCOM-BC test results for each taxon,
    including W-statistic, p-value, adjusted p-value, and whether or not the result is statistically
    significant. See `skbio.stats.composition.ancombc` for details.
    - `reference_group`: Group name used as the reference. All other groups are compared to this
    reference group.
    - `adjustment_method`: method used to adjust p-values to control the false discovery rate
    - `alpha`: p-value threshold used to determine statistical significance
    - `sample_size`: number of samples used in the test after filtering
    - `group_by_variable`: name of the variable used to group samples by
    - `group_sizes`: dict mapping group name to sample size in each group
    """

    main_results: pd.DataFrame
    global_results: pd.DataFrame | None = None
    reference_group: str
    adjustment_method: AdjustmentMethod


class StatsMixin(DistanceMixin, BaseSampleCollection):
    def alpha_diversity_stats(
        self,
        *,
        group_by: str | tuple[str, ...] | list[str],
        paired_by: Optional[str | tuple[str, ...] | list[str]] = None,
        metric: Metric = Metric.Auto,
        test: AlphaDiversityStatsTest = AlphaDiversityStatsTest.Auto,
        diversity_metric: AlphaDiversityMetric = AlphaDiversityMetric.Shannon,
        rank: Rank = Rank.Auto,
        alpha: float = 0.05,
        require_classification_version_match: bool = True,
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

        test : :class:`~onecodex.lib.enums.AlphaDiversityStatsTest`, optional
            Stats test to perform. If `'auto'`, `'mannwhitneyu'` will be chosen if there are two
            groups of unpaired data. `'wilcoxon'` will be chosen if there are two groups and
            `paired_by` is specified. `'kruskal'` will be chosen if there are more than 2 groups.

        rank : :class:`~onecodex.lib.enums.Rank`, optional
            Analysis will be restricted to abundances of taxa at the specified level.
            See :class:`~onecodex.lib.enums.Rank` for details.

        metric: :class:`~onecodex.lib.enums.Metric`, optional
            The taxonomic abundance metric to use. See :class:`~onecodex.lib.enums.Metric`
            for definitions.

        diversity_metric : :class:`~onecodex.lib.enums.AlphaDiversityMetric`
            Function to use when calculating the distance between two samples.

        alpha : float, optional
            Threshold to determine statistical significance when `test="kruskal"`
            (e.g. p < `alpha`). Must be between 0 and 1 (exclusive). If the Kruskal-Wallis p-value
            is significant and there are more than two groups, a posthoc Dunn test is performed.

        require_classification_version_match : bool, optional
            If ``True``, require the same primary classification job ID across all samples included
            in the test.

        Returns
        -------
        :class:`~onecodex.stats.AlphaDiversityStatsResults`

        See Also
        --------
        scipy.stats.wilcoxon
        scipy.stats.mannwhitneyu
        scipy.stats.kruskal
        scikit_posthocs.posthoc_dunn

        """
        self._assert_valid_alpha(alpha)

        metric, rank = self._parse_classification_config_args(metric=metric, rank=rank)

        group_by = self._tuplize(group_by)
        metadata_fields = [group_by]
        if paired_by is not None:
            paired_by = self._tuplize(paired_by)
            metadata_fields.append(paired_by)

        # Munge the metadata first in case there's any errors
        metadata_results = self._metadata_fetch(
            metadata_fields,
            results_df=self.to_classification_df(rank=rank, metric=metric),
            coerce_missing_composite_fields=False,
        )

        df = metadata_results.df
        magic_fields = metadata_results.renamed_fields
        group_by_column_name = magic_fields[group_by]
        paired_by_column_name = magic_fields.get(paired_by, None)

        # Compute alpha diversity
        df.loc[:, diversity_metric] = self.alpha_diversity(
            diversity_metric=diversity_metric, metric=metric, rank=rank
        )

        # Drop NaN alpha diversity values
        prev_count = len(df)
        df.dropna(subset=[diversity_metric], inplace=True)
        num_dropped = prev_count - len(df)
        if num_dropped > 0:
            warnings.warn(
                f"{num_dropped} sample{'s were' if num_dropped > 1 else ' was'} excluded from the "
                f"test as no {diversity_metric} value was calculated for "
                f"{'them' if num_dropped > 1 else 'it'}.",
                StatsWarning,
            )

        # Drop samples with missing group_by data
        df = self._drop_missing_data(df, group_by_column_name, "group_by")

        if paired_by_column_name:
            # Drop samples with missing paired_by data
            df = self._drop_missing_data(df, paired_by_column_name, "paired_by")

        # Drop samples missing abundance data
        if metric.is_abundance_metric:
            df = self._drop_classifications_without_abundances(df)

        # Drop groups of size < 2
        # NOTE: it's important to do this filtering *after* the other filters in case those filters
        # change the group sizes.
        df = self._drop_group_sizes_smaller_than(df, group_by_column_name, 2)

        self._assert_min_num_groups(df, group_by_column_name, 2)

        if require_classification_version_match:
            self._assert_single_database_version(df)

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
            return self._wilcoxon(
                df, diversity_metric, group_by_column_name, paired_by_column_name, alpha
            )
        elif test == AlphaDiversityStatsTest.Mannwhitneyu:
            return self._mannwhitneyu(df, diversity_metric, group_by_column_name, alpha)
        elif test == AlphaDiversityStatsTest.Kruskal:
            return self._kruskal(df, diversity_metric, group_by_column_name, alpha)
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

    def _assert_single_database_version(self, df: pd.DataFrame):
        classification_ids = set(df.index)
        job_ids = {c.job.id for c in self._classifications if c.id in classification_ids}
        if len(job_ids) > 1:
            raise StatsException(
                "To prevent the possibility of false results due to differences in database "
                "versions, running statistical tests is disabled for collections of samples "
                "with different database versions as their primary classifications."
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

    def _get_group_sizes(self, df: pd.DataFrame, group_by_column_name: str) -> dict[str, int]:
        return {
            name: len(group_df)
            for name, group_df in df.groupby(group_by_column_name, observed=True)
        }

    def _get_group_sizes_and_alpha_diversity_values(
        self, df: pd.DataFrame, metric: str, group_by_column_name: str
    ) -> tuple[dict[str, int], list[pd.Series]]:
        group_sizes = {}
        group_alpha_values = []
        for group_name, group_df in df.groupby(group_by_column_name):
            group_sizes[group_name] = len(group_df)
            group_alpha_values.append(group_df[metric])

        return group_sizes, group_alpha_values

    def _wilcoxon(
        self,
        df: pd.DataFrame,
        diversity_metric: AlphaDiversityMetric,
        group_by_column_name: str,
        paired_by_column_name: Optional[str],
        alpha: float,
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

        result = wilcoxon(group1_df[diversity_metric], group2_df[diversity_metric])

        return AlphaDiversityStatsResults(
            test=AlphaDiversityStatsTest.Wilcoxon,
            statistic=result.statistic,
            pvalue=result.pvalue,
            alpha=alpha,
            sample_size=len(df),
            group_by_variable=group_by_column_name,
            group_sizes={group1_name: len(group1_df), group2_name: len(group2_df)},
            paired_by_variable=paired_by_column_name,
        )

    def _mannwhitneyu(
        self,
        df: pd.DataFrame,
        diversity_metric: AlphaDiversityMetric,
        group_by_column_name: str,
        alpha: float,
    ) -> AlphaDiversityStatsResults:
        from scipy.stats import mannwhitneyu

        self._assert_exact_num_groups(df, group_by_column_name, 2)

        group_sizes, group_alpha_values = self._get_group_sizes_and_alpha_diversity_values(
            df, diversity_metric, group_by_column_name
        )
        result = mannwhitneyu(*group_alpha_values)

        return AlphaDiversityStatsResults(
            test=AlphaDiversityStatsTest.Mannwhitneyu,
            statistic=result.statistic,
            pvalue=result.pvalue,
            alpha=alpha,
            sample_size=len(df),
            group_by_variable=group_by_column_name,
            group_sizes=group_sizes,
        )

    def _kruskal(
        self,
        df: pd.DataFrame,
        diversity_metric: AlphaDiversityMetric,
        group_by_column_name: str,
        alpha: float,
    ) -> AlphaDiversityStatsResults:
        from scikit_posthocs import posthoc_dunn
        from scipy.stats import kruskal

        group_sizes, group_alpha_values = self._get_group_sizes_and_alpha_diversity_values(
            df, diversity_metric, group_by_column_name
        )
        result = kruskal(*group_alpha_values)

        posthoc = None
        if result.pvalue < alpha and len(group_sizes) > 2:
            posthoc = PosthocResults(
                test=PosthocStatsTest.Dunn,
                adjustment_method=AdjustmentMethod.BenjaminiHochberg,
                adjusted_pvalues=posthoc_dunn(
                    df, val_col=diversity_metric, group_col=group_by_column_name, p_adjust="fdr_bh"
                ),
            )

        return AlphaDiversityStatsResults(
            test=AlphaDiversityStatsTest.Kruskal,
            statistic=result.statistic,
            pvalue=result.pvalue,
            alpha=alpha,
            sample_size=len(df),
            group_by_variable=group_by_column_name,
            group_sizes=group_sizes,
            posthoc=posthoc,
        )

    def beta_diversity_stats(
        self,
        *,
        group_by: str | tuple[str, ...] | list[str],
        metric: Metric = Metric.Auto,
        diversity_metric: BetaDiversityMetric = BetaDiversityMetric.BrayCurtis,
        rank: Rank = Rank.Auto,
        alpha: float = 0.05,
        num_permutations: int = 999,
        require_classification_version_match: bool = True,
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

        metric: :class:`~onecodex.lib.enums.Metric`, optional
            The taxonomic abundance metric to use. See :class:`~onecodex.lib.enums.Metric`
            for definitions.

        diversity_metric : :class:`~onecodex.lib.enums.BetaDiversityMetric`
            Function to use when calculating the distance between two samples.

        rank : :class:`~onecodex.lib.enums.Rank`, optional
            Analysis will be restricted to abundances of taxa at the specified level.
            See :class:`~onecodex.lib.enums.Rank` for details.

        alpha : float, optional
            Threshold to determine statistical significance (e.g. p < `alpha`). Must be between 0
            and 1 (exclusive). If the p-value is significant and there are more than two groups,
            posthoc pairwise PERMANOVA tests are performed.

        num_permutations : int, optional
            Number of permutations to use when computing the p-value.

        require_classification_version_match : bool, optional
            If ``True``, require the same primary classification job ID across all samples included
            in the test.

        Returns
        -------
        :class:`~onecodex.stats.BetaDiversityStatsResults`

        See Also
        --------
        skbio.stats.distance.permanova
        scipy.stats.false_discovery_control

        """
        self._assert_valid_alpha(alpha)
        metric, rank = self._parse_classification_config_args(metric=metric, rank=rank)

        # Munge the metadata first in case there's any errors
        group_by = self._tuplize(group_by)
        metadata_results = self._metadata_fetch(
            [group_by],
            results_df=self.to_classification_df(rank=rank, metric=metric),
            coerce_missing_composite_fields=False,
        )
        df = metadata_results.df
        magic_fields = metadata_results.renamed_fields
        group_by_column_name = magic_fields[group_by]

        # Drop samples with missing group_by data
        df = self._drop_missing_data(df, group_by_column_name, "group_by")

        # Drop samples missing abundance data
        if metric.is_abundance_metric:
            df = self._drop_classifications_without_abundances(df)

        # Drop groups of size < 2
        # NOTE: it's important to do this filtering *after* the other filters in case those filters
        # change the group sizes.
        df = self._drop_group_sizes_smaller_than(df, group_by_column_name, 2)

        self._assert_min_num_groups(df, group_by_column_name, 2)

        if require_classification_version_match:
            self._assert_single_database_version(df)

        # Compute beta diversity distance matrix and filter it to match the remaining IDs in the
        # metadata dataframe
        dm = self.beta_diversity(diversity_metric=diversity_metric, rank=rank).filter(
            df.index, strict=True
        )

        return self._permanova(
            dm=dm,
            df=df,
            group_by_column_name=group_by_column_name,
            alpha=alpha,
            num_permutations=num_permutations,
        )

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

        group_sizes = self._get_group_sizes(df, group_by_column_name)
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
            group_names = sorted(group_sizes)
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
                test=PosthocStatsTest.PairwisePermanova,
                adjustment_method=AdjustmentMethod.BenjaminiHochberg,
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
            alpha=alpha,
            num_permutations=result["number of permutations"],
            sample_size=result["sample size"],
            group_by_variable=group_by_column_name,
            group_sizes=group_sizes,
            posthoc=posthoc,
        )

    def _ancombc(
        self,
        *,
        group_by: str | tuple[str, ...] | list[str],
        reference_group: str | None = None,
        metric: Metric = Metric.Auto,
        rank: Rank = Rank.Auto,
        alpha: float = 0.05,
        include_global_test: bool = False,
        require_classification_version_match: bool = True,
    ) -> AncombcResults:
        """Perform a test for differentially abundant taxa using ANCOM-BC.

        Readcounts are normalized and zeros are replaced using multiplicative replacement. Results
        indicate taxa that are significantly different in abundance across the grouping variable of
        interest, as well as each taxon's log2-fold change (with standard error), W-statistic,
        p-value, and adjusted p-value.

        This method is **experimental** and breaking API changes may occur.

        Parameters
        ----------
        group_by : str or tuple of str or list of str
            Metadata variable to group samples by. At least two groups are required. If `group_by`
            is a tuple or list, field values are joined with an underscore character ("_").

        reference_group : str, optional
            Metadata variable group name to use as the reference. All other groups will be compared
            to this reference group. If not specified, the first alphabetically sorted group name
            will be used as the reference.

        metric: :class:`~onecodex.lib.enums.Metric`, optional
            The taxonomic abundance metric to use. See :class:`~onecodex.lib.enums.Metric`
            for definitions.

        rank : :class:`~onecodex.lib.enums.Rank`, optional
            Analysis will be restricted to abundances of taxa at the specified level.
            See :class:`~onecodex.lib.enums.Rank` for details.

        alpha : float, optional
            Threshold to determine statistical significance (e.g. p < `alpha`). Must be between 0
            and 1 (exclusive).

        include_global_test : bool, optional
            If ``True``, include global ANCOM-BC results in addition to the main results. A global
            test can only be run if there are at least 3 groups.

        require_classification_version_match : bool, optional
            If ``True``, require the same primary classification job ID across all samples included
            in the test.

        Returns
        -------
        :class:`~onecodex.stats.AncombcResults`

        See Also
        --------
        skbio.stats.composition.ancombc
        skbio.stats.composition.multi_replace

        """
        import pandas as pd
        from skbio.stats.composition import ancombc, multi_replace

        self._assert_valid_alpha(alpha)
        metric, rank = self._parse_classification_config_args(metric=metric, rank=rank)

        # Munge the metadata first in case there's any errors
        group_by = self._tuplize(group_by)
        metadata_results = self._metadata_fetch(
            [group_by],
            results_df=self.to_classification_df(rank=rank, metric=metric),
            coerce_missing_composite_fields=False,
        )
        metadata_df = metadata_results.df
        group_by_column_name = metadata_results.renamed_fields[group_by]

        # Drop samples with missing group_by data
        metadata_df = self._drop_missing_data(metadata_df, group_by_column_name, "group_by")

        # Drop samples missing abundance data
        if metric.is_abundance_metric:
            metadata_df = self._drop_classifications_without_abundances(metadata_df)

        taxa_df = self.to_classification_df(rank=rank, metric=metric)
        taxa_df = self._rename_tax_ids(taxa_df)

        # `multi_replace()` will error on rows that are all zeros
        metadata_df = self._drop_samples_with_zero_abundance(metadata_df, taxa_df)

        # Drop groups of size < 2
        # NOTE: it's important to do this filtering *after* the other filters in case those filters
        # change the group sizes
        metadata_df = self._drop_group_sizes_smaller_than(metadata_df, group_by_column_name, 2)

        self._assert_min_num_groups(metadata_df, group_by_column_name, 2)

        if require_classification_version_match:
            self._assert_single_database_version(metadata_df)

        group_names = sorted(metadata_df[group_by_column_name].unique())
        if reference_group:
            if reference_group not in group_names:
                raise StatsException(
                    f"`reference_group` must be a group name defined by `group_by`. Valid group "
                    f"names after filtering: {', '.join(group_names)}"
                )
        else:
            reference_group = metadata_df[group_by_column_name].min()

        # Reference group must be the first group
        group_names.remove(reference_group)
        group_names = [reference_group] + group_names
        metadata_df[group_by_column_name] = pd.Categorical(
            metadata_df[group_by_column_name], categories=group_names, ordered=True
        )

        # ANCOM-BC can only run a global test if there's at least 3 groups
        if include_global_test:
            self._assert_min_num_groups_for_global_test(metadata_df, group_by_column_name)

        # Filter the taxa dataframe to match the remaining IDs in the metadata dataframe
        taxa_df = taxa_df.loc[metadata_df.index]

        # Need to use pseudocounts or `multi_replace()` because ANCOM-BC can't handle zeros
        multi_replace_df = pd.DataFrame(
            multi_replace(taxa_df), columns=taxa_df.columns, index=taxa_df.index
        )

        # Use a safe internal column name for the patsy formula so that it works as both a
        # valid formula term and a metadata column lookup (skbio's `grouping` parameter is
        # used for both). Using patsy's quoting, e.g. `Q('foo bar')`, doesn't work with the
        # `grouping` parameter.
        ancombc_col = "__ancombc_group__"
        ancombc_metadata = metadata_df.rename(columns={group_by_column_name: ancombc_col})
        ancombc_results = ancombc(
            multi_replace_df,
            ancombc_metadata,
            formula=ancombc_col,
            grouping=ancombc_col if include_global_test else None,
            alpha=alpha,
            p_adjust="fdr_bh",
        )
        main_results = ancombc_results
        global_results = None
        if include_global_test:
            main_results, global_results = ancombc_results
            global_results = global_results.rename_axis(index={"FeatureID": "Taxon"})

        # Rename the __ancombc_group__ value in the Covariate column back to the original column
        # name, and reformat the value to be easier to read
        main_results = main_results.rename_axis(index={"FeatureID": "Taxon"}).rename(
            index=lambda covariate: self._rename_ancombc_covariate(
                covariate, ancombc_col, group_by_column_name, reference_group
            ),
            level="Covariate",
        )

        return AncombcResults(
            main_results=main_results,
            global_results=global_results,
            reference_group=reference_group,
            alpha=alpha,
            adjustment_method=AdjustmentMethod.BenjaminiHochberg,
            sample_size=len(metadata_df),
            group_by_variable=group_by_column_name,
            group_sizes=self._get_group_sizes(metadata_df, group_by_column_name),
        )

    def _rename_tax_ids(self, taxa_df: pd.DataFrame) -> pd.DataFrame:
        return taxa_df.rename(
            columns=lambda tax_id: f"{self.taxonomy['name'][tax_id]} ({tax_id})"
            if tax_id in self.taxonomy["name"]
            else tax_id
        )

    def _drop_samples_with_zero_abundance(
        self, metadata_df: pd.DataFrame, taxa_df: pd.DataFrame
    ) -> pd.DataFrame:
        prev_count = len(metadata_df)
        metadata_df = metadata_df.loc[(taxa_df != 0).any(axis=1)]
        num_dropped = prev_count - len(metadata_df)
        if num_dropped > 0:
            warnings.warn(
                f"{num_dropped} sample{'s were' if num_dropped > 1 else ' was'} excluded from the "
                f"test as {'they' if num_dropped > 1 else 'it'} had zero abundance across all "
                f"taxa.",
                StatsWarning,
            )
        return metadata_df

    def _assert_min_num_groups_for_global_test(
        self, metadata_df: pd.DataFrame, group_by_column_name: str
    ):
        try:
            self._assert_min_num_groups(metadata_df, group_by_column_name, 3)
        except StatsException:
            # Raise a more specific error message
            raise StatsException(
                "`include_global_test=True` may only be used if there are at least 3 groups after "
                "filtering."
            )

    def _rename_ancombc_covariate(
        self,
        covariate: str,
        ancombc_col: str,
        group_by_column_name: str,
        reference_group: str,
    ) -> str:
        # Format is "__ancombc_group__[T.<group_name>]". Reformat to
        # "<group_by_column_name>: <group_name> vs <reference_group> (reference)"
        m = re.fullmatch(rf"{re.escape(ancombc_col)}\[T\.(.+)\]", covariate)
        if m:
            return f"{group_by_column_name}: {m.group(1)} vs {reference_group} (reference)"
        return covariate
