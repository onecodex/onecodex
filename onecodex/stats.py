from __future__ import annotations
from dataclasses import dataclass
import warnings
from typing import TYPE_CHECKING, Optional

from onecodex.exceptions import StatsException, StatsWarning, PlottingException
from onecodex.lib.enums import Rank
from onecodex.viz._primitives import prepare_props

if TYPE_CHECKING:
    import altair as alt
    import pandas as pd


@dataclass(frozen=True)
class DifferentialAbundanceResults:
    # ANCOM results
    ancom_df: pd.DataFrame
    percentile_df: pd.DataFrame

    # Intermediate data needed for plotting
    _taxa_df: pd.DataFrame
    _metadata_df: pd.DataFrame
    _group_by_column_name: str
    _metric_name: str

    def plot(
        self,
        title: Optional[str] = None,
        width: Optional[int | str] = None,
        height: Optional[int | str] = None,
        return_chart: bool = False,
    ) -> Optional[alt.FacetChart]:
        """Plot differential abundance results as boxplots faceted by taxa.

        A boxplot of abundances per group is generated for each significantly different taxon.

        Parameters
        ----------
        title : str, optional
            Text label at the top of the chart.
        width : int or str, optional
            The width of the chart. Can be `"container"` for responsive sizing.
        height : int or str, optional
            The height of the chart. Can be `"container"` for responsive sizing.
        return_chart : bool, optional
            When `True`, return an `altair.FacetChart` object instead of displaying the resulting
            chart in the current notebook.

        Returns
        -------
        alt.FacetChart or None
            If `return_chart` is `True`, the Altair chart is returned, otherwise `None` is returned
            and the chart is displayed in the current notebook.

        """
        import altair as alt
        import pandas as pd

        reject_column = self.ancom_df["Reject null hypothesis"]
        if not reject_column.any():
            raise PlottingException(
                "No significantly different taxa were detected in ANCOM results."
            )

        # Filter to taxa that are significantly different
        filtered_taxa_df = self._taxa_df[reject_column[reject_column].index].copy()

        # Add grouping metadata variable
        filtered_taxa_df.loc[:, self._group_by_column_name] = self._metadata_df[
            self._group_by_column_name
        ]

        # - Facet based on the different taxa that were found to be significantly different
        # - Color by the grouping metadata variable
        plotting_df = pd.melt(
            filtered_taxa_df,
            id_vars=self._group_by_column_name,
            var_name="Significantly Different Taxa",
            value_name=self._metric_name,
        )
        chart = (
            alt.Chart(plotting_df)
            .mark_boxplot(median={"stroke": "black"})
            .encode(
                x=alt.X(self._group_by_column_name),
                y=alt.Y(self._metric_name),
                color=self._group_by_column_name,
            )
            .facet("Significantly Different Taxa")
            .resolve_scale(y="independent")
        )
        chart = chart.properties(**prepare_props(title=title, width=width, height=height))

        return chart if return_chart else chart.display()


class StatsMixin:
    def differential_abundance(
        self,
        *,
        group_by: str | tuple[str, ...],
        rank: Rank = Rank.Auto,
        alpha: float = 0.05,
    ) -> DifferentialAbundanceResults:
        """
        Perform a test for differentially abundant taxa using ANCOM.

        Readcounts are normalized and zeros are replaced using multiplicative replacement. Results
        indicate taxa that are significantly different in abundance across the grouping variable of
        interest, as well as each taxon's W-statistic.

        Parameters
        ----------
        group_by : str or tuple of str
            Metadata variable to group samples by. At least two groups are required. If `group_by`
            is a tuple, field values are joined with an underscore character ("_").
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.
        alpha : float, optional
            Threshold to determine statistical significance (e.g. p < `alpha`). Must be between 0
            and 1 (exclusive).

        Returns
        -------
        DifferentialAbundanceResults
            A dataclass with these attributes:

            - `ancom_df`: `pd.DataFrame` containing ANCOM test results: taxa, their W-statistics,
            and whether or not they are significantly different. See `skbio.stats.composition.ancom`
            for details.

            - `percentile_df`: `pd.DataFrame` containing percentile abundances of taxa in each
            group. See `skbio.stats.composition.ancom` for details.

            There is also a `plot()` method for generating faceted boxplots.

        See Also
        --------
        skbio.stats.composition.ancom
        skbio.stats.composition.multiplicative_replacement

        """
        import pandas as pd
        from skbio.stats.composition import ancom, multiplicative_replacement

        if not isinstance(group_by, tuple):
            group_by = (group_by,)

        # Munge the metadata first in case there's any errors
        metadata_df, magic_fields = self._metadata_fetch([group_by])

        # Normalize because `multiplicative_replacement()` requires compositions
        taxa_df = self.to_df(rank=rank, normalize=True)
        taxa_df = taxa_df.rename(
            columns=lambda tax_id: f"{self.taxonomy['name'][tax_id]} ({tax_id})"
            if tax_id in self.taxonomy["name"]
            else tax_id
        )

        # `multiplicative_replacement()` will error on rows that are all zeros
        prev_index = taxa_df.index
        taxa_df = taxa_df.loc[(taxa_df != 0).any(axis=1)]
        num_dropped = len(prev_index) - len(taxa_df.index)
        if num_dropped > 0:
            warnings.warn(
                f"{num_dropped} sample{'s' if num_dropped > 1 else ''} with zero abundance across "
                f"all taxa {'were' if num_dropped > 1 else 'was'} ignored.",
                StatsWarning,
            )

        # Drop any samples from the metadata df that are not in the taxa df
        metadata_df = metadata_df.filter(items=taxa_df.index, axis="index")

        group_by_column_name = magic_fields[group_by]
        num_groups = metadata_df[group_by_column_name].nunique(dropna=False)
        if num_groups < 2:
            raise StatsException("`group_by` must have at least two groups to test between.")
        if num_groups == len(metadata_df):
            raise StatsException("Each group defined by `group_by` contains only one value.")

        # Need to use pseudocounts or `multiplicative_replacement()` because ANCOM can't handle
        # zeros
        multiplicative_replacement_df = pd.DataFrame(
            multiplicative_replacement(taxa_df), columns=taxa_df.columns, index=taxa_df.index
        )
        ancom_df, percentile_df = ancom(
            multiplicative_replacement_df,
            metadata_df[group_by_column_name],
            alpha=alpha,
            multiple_comparisons_correction="holm-bonferroni",
        )

        return DifferentialAbundanceResults(
            ancom_df=ancom_df,
            percentile_df=percentile_df,
            _taxa_df=taxa_df,
            _metadata_df=metadata_df,
            _group_by_column_name=group_by_column_name,
            _metric_name=self.metric,
        )
