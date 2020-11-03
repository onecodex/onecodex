from onecodex.lib.enums import Rank, Linkage
from onecodex.exceptions import OneCodexException, PlottingException
from onecodex.viz._primitives import prepare_props, sort_helper, get_base_classification_url


class VizHeatmapMixin(object):
    def plot_heatmap(
        self,
        rank=Rank.Auto,
        normalize="auto",
        top_n="auto",
        threshold="auto",
        title=None,
        xlabel=None,
        ylabel=None,
        tooltip=None,
        return_chart=False,
        linkage=Linkage.Average,
        haxis=None,
        metric="euclidean",
        legend="auto",
        label=None,
        sort_x=None,
        sort_y=None,
        width=None,
        height=None,
    ):
        """Plot heatmap of taxa abundance/count data for several samples.

        Parameters
        ----------
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.
        normalize : 'auto' or `bool`, optional
            Convert read counts to relative abundances such that each sample sums to 1.0. Setting
            'auto' will choose automatically based on the data.
        return_chart : `bool`, optional
            When True, return an `altair.Chart` object instead of displaying the resulting plot in
            the current notebook.
        haxis : `string`, optional
            The metadata field (or tuple containing multiple categorical fields) used to group
            samples together. Each group of samples will be clustered independently.
        metric : {'euclidean', 'braycurtis', 'cityblock', 'manhattan', 'jaccard', 'unifrac', 'unweighted_unifrac', 'aitchison'}, optional
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

        Examples
        --------
        Plot a heatmap of the relative abundances of the top 10 most abundant families.

        >>> plot_heatmap(rank='family', top_n=10)
        """
        # Deferred imports
        import altair as alt
        import pandas as pd

        if rank is None:
            raise OneCodexException("Please specify a rank or 'auto' to choose automatically")

        if not (threshold or top_n):
            raise OneCodexException("Please specify at least one of: threshold, top_n")

        if len(self._results) < 2:
            raise PlottingException(
                "There are too few samples for heatmap plots after filtering. Please select 2 or "
                "more samples to plot."
            )

        if top_n == "auto" and threshold == "auto":
            top_n = 10
            threshold = None
        elif top_n == "auto" and threshold != "auto":
            top_n = None
        elif top_n != "auto" and threshold == "auto":
            threshold = None

        df = self.to_df(
            rank=rank, normalize=normalize, top_n=top_n, threshold=threshold, table_format="long"
        )

        if len(df["tax_id"].unique()) < 2:
            raise PlottingException(
                "There are too few taxa for heatmap clustering after filtering. Please select a "
                "rank or threshold that includes at least 2 taxa."
            )

        if legend == "auto":
            legend = df.ocx.metric

        if tooltip:
            if not isinstance(tooltip, list):
                tooltip = [tooltip]
        else:
            tooltip = []

        if haxis:
            tooltip.append(haxis)

        tooltip.insert(0, "Label")

        magic_metadata, magic_fields = self._metadata_fetch(tooltip, label=label)

        # add columns for prettier display
        df["Label"] = magic_metadata["Label"][df["classification_id"]].tolist()
        df["tax_name"] = ["{} ({})".format(self.taxonomy["name"][t], t) for t in df["tax_id"]]

        # and for metadata
        for f in tooltip:
            df[magic_fields[f]] = magic_metadata[magic_fields[f]][df["classification_id"]].tolist()

        # if we've already been normalized, we must cluster samples by euclidean distance. beta
        # diversity measures won't work with normalized distances.
        if self._guess_normalized():
            if metric != "euclidean":
                raise OneCodexException(
                    "Results are normalized. Please re-run with metric=euclidean"
                )

            df_sample_cluster = self.to_df(
                rank=rank, normalize=normalize, top_n=top_n, threshold=threshold
            )
            df_taxa_cluster = df_sample_cluster
        else:
            df_sample_cluster = self.to_df(
                rank=rank, normalize=False, top_n=top_n, threshold=threshold
            )

            df_taxa_cluster = self.to_df(
                rank=rank, normalize=normalize, top_n=top_n, threshold=threshold
            )

        # applying clustering to determine order of taxa, or use custom sorting function if given
        if sort_y is None:
            taxa_cluster = df_taxa_cluster.ocx._cluster_by_taxa(linkage=linkage)
            taxa_cluster = taxa_cluster["labels_in_order"]
        else:
            taxa_cluster = sort_helper(sort_y, df["tax_name"])

        if sort_x is None:
            if haxis is None:
                # cluster samples only once
                sample_cluster = df_sample_cluster.ocx._cluster_by_sample(
                    rank=rank, metric=metric, linkage=linkage
                )
                labels_in_order = magic_metadata["Label"][sample_cluster["ids_in_order"]].tolist()
            else:
                if not (
                    pd.api.types.is_bool_dtype(df[magic_fields[haxis]])
                    or pd.api.types.is_categorical_dtype(df[magic_fields[haxis]])  # noqa
                    or pd.api.types.is_object_dtype(df[magic_fields[haxis]])  # noqa
                ):  # noqa
                    raise OneCodexException(
                        "Metadata field on horizontal axis can not be numerical"
                    )

                labels_in_order = []
                df_sample_cluster[haxis] = self.metadata[haxis]

                for group, group_df in df_sample_cluster.groupby(haxis):

                    if group_df.shape[0] <= 3:
                        # we can't cluster
                        labels_in_order.extend(
                            sorted(magic_metadata["Label"][group_df.index].tolist())
                        )
                        continue

                    sample_cluster = group_df.drop(columns=[haxis]).ocx._cluster_by_sample(
                        rank=rank, metric=metric, linkage=linkage
                    )
                    labels_in_order.extend(
                        magic_metadata["Label"][sample_cluster["ids_in_order"]].tolist()
                    )
        else:
            labels_in_order = sort_helper(sort_x, magic_metadata["Label"].tolist())

        # should ultimately be Label, tax_name, readcount_w_children, then custom fields
        tooltip_for_altair = [magic_fields[f] for f in tooltip]
        tooltip_for_altair.insert(1, "tax_name")
        tooltip_for_altair.insert(2, "{}:Q".format(df.ocx.metric))

        alt_kwargs = dict(
            x=alt.X("Label:N", axis=alt.Axis(title=xlabel), sort=labels_in_order),
            y=alt.Y("tax_name:N", axis=alt.Axis(title=ylabel), sort=taxa_cluster),
            color=alt.Color("{}:Q".format(df.ocx.metric), legend=alt.Legend(title=legend)),
            tooltip=tooltip_for_altair,
            href="url:N",
            url=get_base_classification_url() + alt.datum.classification_id,
        )

        if haxis:
            alt_kwargs["column"] = alt.Column(
                haxis, header=alt.Header(titleOrient="bottom", labelOrient="bottom")
            )

        chart = (
            alt.Chart(df)
            .transform_calculate(url=alt_kwargs.pop("url"))
            .mark_rect()
            .encode(**alt_kwargs)
        )

        col_count = len(labels_in_order)
        row_count = len(taxa_cluster)

        chart = chart.properties(
            **prepare_props(
                title=title, height=(height or 15 * row_count), width=(width or 15 * col_count)
            )
        )

        if haxis:
            chart = chart.resolve_scale(x="independent")

        if return_chart:
            return chart
        else:
            chart.interactive().display()
