import warnings
from typing import Union

from onecodex.exceptions import OneCodexException, PlottingException, PlottingWarning
from onecodex.lib.enums import (
    BetaDiversityMetric,
    Link,
    Linkage,
    Metric,
    Rank,
)
from onecodex.models.base_sample_collection import BaseSampleCollection
from onecodex.viz._primitives import (
    escape_chart_fields,
    get_classification_url,
    get_ncbi_taxonomy_browser_url,
    prepare_props,
    sort_helper,
)


class VizHeatmapMixin(BaseSampleCollection):
    def plot_heatmap(
        self,
        rank: Union[Rank, str] = Rank.Auto,
        metric: Union[Metric, str] = Metric.Auto,
        top_n="auto",
        threshold="auto",
        title=None,
        xlabel=None,
        ylabel=None,
        tooltip=None,
        return_chart=False,
        linkage=Linkage.Average,
        haxis=None,
        diversity_metric: BetaDiversityMetric = BetaDiversityMetric.Euclidean,
        legend="auto",
        label=None,
        sort_x=None,
        sort_y=None,
        width=None,
        height=None,
        link=Link.Ocx,
        match_taxonomy=True,
    ):
        """Plot heatmap of taxa abundance/count data for several samples.

        Parameters
        ----------
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.
        metric: Metric
            the taxonomic abundance metric to use. See onecodex.lib.enums.Metric for definitions
        return_chart : `bool`, optional
            When True, return an `altair.Chart` object instead of displaying the resulting plot in
            the current notebook.
        haxis : `string`, optional
            The metadata field (or tuple containing multiple categorical fields) used to group
            samples together. Each group of samples will be clustered independently.
        diversity_metric: BetaDiversityMetric
            Function to use when calculating the distance between two samples.
            Note that 'cityblock' and 'manhattan' are equivalent metrics.
        linkage : {'average', 'single', 'complete', 'weighted', 'centroid', 'median'}
            The type of linkage to use when clustering axes.
        top_n : `int`, optional
            Display the top N most abundant taxa in the entire cohort of samples.
        threshold : `float`
            Display only taxa that are more abundant that this threshold in one or more samples.
        title : `string`, optional
            Text label at the top of the plot.
        xlabel : `string`, optional
            Text label along the horizontal axis.
        ylabel : `string`, optional
            Text label along the vertical axis.
        tooltip : `string` or `list`, optional
            A string or list containing strings representing metadata fields. When a point in the
            plot is hovered over, the value of the metadata associated with that sample will be
            displayed in a modal.
        legend: `string`, optional
            Title for color scale. Defaults to the field used to generate the plot, e.g.
            readcount_w_children or abundance.
        label : `string` or `callable`, optional
            A metadata field (or function) used to label each analysis. If passing a function, a
            dict containing the metadata for each analysis is passed as the first and only
            positional argument. The callable function must return a string.
        sort_x : `list` or `callable`, optional
            Either a list of sorted labels or a function that will be called with a list of x-axis labels
            as the only argument, and must return the same list in a user-specified order.
        sort_y : `list` or `callable`, optional
            Either a list of sorted labels or a function that will be called with a list of y-axis labels
            as the only argument, and must return the same list in a user-specified order.
        link : {'ocx', 'ncbi'}, optional
            If `link` is 'ocx', clicking a sample will open its classification results in the One
            Codex app. If `link` is 'ncbi', clicking a taxon will open the NCBI taxonomy browser.
        match_taxonomy : `bool`, default=True
            Whether or not to consider taxonomic names when looking for metadata fields mapped to
            plot attributes including `tooltip`, `label`

        Examples
        --------
        Plot a heatmap of the relative abundances of the top 10 most abundant families.

        >>> samples.plot_heatmap(rank='family', top_n=10)
        """
        # Deferred imports
        import altair as alt
        import numpy as np
        import pandas as pd

        metric, rank = self._parse_classification_config_args(metric=metric, rank=rank)

        if not (threshold or top_n):
            raise OneCodexException("Please specify at least one of: threshold, top_n")

        if len(self._classifications) < 2:
            raise PlottingException(
                "There are too few samples for heatmap plots after filtering. Please select 2 or "
                "more samples to plot."
            )

        if metric.is_abundance_metric and (
            len(self._classifications) == len(self._classification_ids_without_abundances)
        ):
            raise PlottingException(
                "Abundances are not calculated for any of the selected samples. Please select a "
                "different metric or a different set of samples to plot."
            )

        if top_n == "auto" and threshold == "auto":
            top_n = 10
            threshold = None
        elif top_n == "auto" and threshold != "auto":
            top_n = None
        elif top_n != "auto" and threshold == "auto":
            threshold = None

        df = self.to_classification_df(
            rank=rank,
            metric=metric,
            top_n=top_n,
            table_format="long",
            threshold=threshold,
            fill_missing=False,
        )

        if len(df["tax_id"].unique()) < 2:
            raise PlottingException(
                "There are too few taxa for heatmap clustering after filtering. Please select a "
                "rank or threshold that includes at least 2 taxa."
            )

        if legend == "auto":
            legend = df.ocx_metric

        if tooltip:
            if not isinstance(tooltip, list):
                tooltip = [tooltip]
        else:
            tooltip = []

        if haxis:
            tooltip.append(haxis)

        tooltip.insert(0, "Label")

        metadata_results = self._metadata_fetch(
            tooltip, results_df=df, label=label, match_taxonomy=match_taxonomy
        )
        magic_metadata = metadata_results.df
        magic_fields = metadata_results.renamed_fields

        magic_metadata.replace(np.nan, "N/A", inplace=True)

        # add columns for prettier display
        df["Label"] = magic_metadata["Label"][df["classification_id"]].tolist()

        df["tax_name"] = ["{} ({})".format(self.taxonomy["name"][t], t) for t in df["tax_id"]]

        # and for metadata
        for f in tooltip:
            df[magic_fields[f]] = magic_metadata[magic_fields[f]][df["classification_id"]].tolist()

        # if we've already been normalized, we must cluster samples by euclidean distance. beta
        # diversity measures won't work with normalized distances.
        if metric.is_normalized and diversity_metric != BetaDiversityMetric.Euclidean:
            warnings.warn(
                f"Only euclidean distances are supported for metric {metric}", PlottingWarning
            )

        # df is long, so create a wide version for clustering
        df_cluster = self.to_df(
            metric=metric, rank=rank, top_n=top_n, threshold=threshold, fill_missing=False
        )

        # applying clustering to determine order of taxa, or use custom sorting function if given
        if sort_y is None:
            taxa_cluster = self._cluster_by_taxa(results_df=df_cluster, linkage=linkage)
            taxa_cluster = taxa_cluster["labels_in_order"]
        else:
            taxa_cluster = sort_helper(sort_y, df["tax_name"])

        if sort_x is None:
            if haxis is None:
                # cluster samples only once
                sample_cluster = self._cluster_by_sample(
                    results_df=df_cluster,
                    diversity_metric=diversity_metric,
                    linkage=linkage,
                )
                labels_in_order = magic_metadata["Label"][sample_cluster["ids_in_order"]].tolist()
            else:
                if not (
                    pd.api.types.is_bool_dtype(df[magic_fields[haxis]])
                    or isinstance(df[magic_fields[haxis]].dtype, pd.CategoricalDtype)
                    or pd.api.types.is_object_dtype(df[magic_fields[haxis]])
                ):
                    raise OneCodexException(
                        "Metadata field on horizontal axis can not be numerical"
                    )

                labels_in_order = []
                df_cluster[haxis] = self.metadata[haxis]

                for group, group_df in df_cluster.groupby(haxis, dropna=False):
                    if group_df.shape[0] <= 3:
                        # we can't cluster
                        labels_in_order.extend(
                            sorted(magic_metadata["Label"][group_df.index].tolist())
                        )
                        continue

                    sample_cluster = self._cluster_by_sample(
                        results_df=group_df.drop(columns=[haxis]),
                        diversity_metric=diversity_metric,
                        linkage=linkage,
                    )
                    labels_in_order.extend(
                        magic_metadata["Label"][sample_cluster["ids_in_order"]].tolist()
                    )
        else:
            labels_in_order = sort_helper(sort_x, magic_metadata["Label"].tolist())

        pretty_metric_name = self._display_name_for_metric(metric=df.ocx_metric)

        # should ultimately be Label, tax_name, readcount_w_children, then custom fields
        tooltip_for_altair = [magic_fields[f] for f in tooltip]
        tooltip_for_altair.insert(1, "tax_name")
        tooltip_for_altair.insert(2, f"{pretty_metric_name}:Q")

        alt_kwargs = dict(
            x=alt.X("Label:N", axis=alt.Axis(title=xlabel), sort=labels_in_order),
            y=alt.Y("tax_name:N", axis=alt.Axis(title=ylabel), sort=taxa_cluster),
            color=alt.Color(
                f"{pretty_metric_name}:Q",
                legend=alt.Legend(title=legend),
            ),
            tooltip=tooltip_for_altair,
        )

        if link == Link.Ocx:
            df["url"] = df["classification_id"].apply(get_classification_url)
            alt_kwargs["href"] = "url:N"
        elif link == Link.Ncbi:
            df["url"] = df["tax_id"].apply(get_ncbi_taxonomy_browser_url)
            alt_kwargs["href"] = "url:N"

        if haxis:
            alt_kwargs["column"] = alt.Column(
                haxis, header=alt.Header(title=haxis, titleOrient="bottom", labelOrient="bottom")
            )

        # Drop all-NaN rows, convert remaining NaNs to 0s, and concat
        # dropped all-NaN rows to the end (only if using Abundance or AbundanceWChildren)
        if metric in {Metric.Abundance, Metric.AbundanceWChildren}:
            end_labels = []
            dropped = []
            for classification_id in self._classification_ids_without_abundances:
                d = df[(df == classification_id).any(axis=1)]
                label = d["Label"].values[0]
                labels_in_order.remove(label)
                end_labels.append(label)
                dropped.append(d)
                index = df[(df == classification_id).any(axis=1)].index
                df = df.drop(index)

            # replace NaN with 0.0 in the not-dropped classifications (classifications that have abundances)
            df = df.replace(np.nan, 0)
            df = pd.concat([df, *dropped])

            labels_in_order = labels_in_order + end_labels

        assert set(df["Label"].values) == set(labels_in_order)

        chart = alt.Chart(df).mark_rect().encode(**alt_kwargs)

        col_count = len(labels_in_order)
        row_count = len(taxa_cluster)

        chart = chart.properties(
            **prepare_props(
                title=title, height=(height or 15 * row_count), width=(width or 15 * col_count)
            )
        )

        if haxis:
            chart = chart.resolve_scale(x="independent")

        escape_chart_fields(chart)

        if return_chart:
            return chart
        else:
            chart.interactive().display()
