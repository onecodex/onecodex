from onecodex.lib.enums import (
    FunctionalAnnotations,
    FunctionalAnnotationsMetric,
)
from onecodex.viz._primitives import prepare_props, sort_helper


class VizFunctionalHeatmapMixin(object):
    def plot_functional_heatmap(
        self,
        num_of_functions=10,
        annotation=FunctionalAnnotations.Go,
        metric=None,
        sort_x=None,
        label=None,
        haxis=None,
        return_chart=False,
        xlabel="",
        title=None,
        width=None,
        height=None,
    ):
        """Plot heatmap of functional profile data.

        Parameters
        ----------
        num_of_functions : `int`, optional
            Display the top N most abundant or covered functions in the entire cohort of samples.
        annotation : `FunctionalAnnotations` or `str`, optional
            {'go', 'eggnog', 'ko', 'ec', 'pfam', 'pathways'}
            Grouping (annotation) sub-database used to group gene families by.
        metric : `FunctionalAnnotationsMetric` or `str`, optional
            {'cpm', 'rpk', 'abundance', 'coverage'}
            Normalization or value to display.
            If annotation is one of 'go', 'eggnog', 'ko', 'ec' or 'pfam', then available metrics include
                'rpk' (read counts normalized by kilobase of gene length), or
                'cpm' (relative copy of gene depth, normalized to a million RPK total).
            If pathways are selected for annotation, then available metrics include
                'abundance' (summed copy numbers of reactions' constituent enzymes)
                'coverage' (probabilistic measure of a complete metabolic pathway,
                    where 1 is high confidence of a complete pathway is covered, and 0 is low confidence)
        sort_x : `list` or `callable`, optional
            Either a list of sorted labels or a function that will be called with a list of x-axis labels
            as the only argument, and must return the same list in a user-specified order.
        label : `str` or `callable`, optional
            A metadata field (or function) used to label each analysis. If passing a function, a
            dict containing the metadata for each analysis is passed as the first and only
            positional argument. The callable function must return a string.
        haxis : `string`, optional
            The metadata field (or tuple containing multiple categorical fields) used to facet
            samples.
        return_chart : `bool`, optional
            When True, return an `altair.Chart` object instead of displaying the resulting plot in
            the current notebook.
        xlabel : `str, optional
            Text label along the horizontal axis.
        title : `str`, optional
            Text label at the top of the plot.
        width : `float` or `str` or `dict`, optional
            Set `altair.Chart.width`.
        height : `float` or `str` or `dict`, optional
            Set `altair.Chart.height`.
        """
        # Deferred imports
        import altair as alt
        import pandas as pd

        # Preparing params
        # ----------------

        annotation = FunctionalAnnotations(annotation)
        if metric is None:
            if annotation == FunctionalAnnotations.Pathways:
                metric = FunctionalAnnotationsMetric.Abundance
            else:
                metric = FunctionalAnnotationsMetric.Cpm
        metric = FunctionalAnnotationsMetric(metric)

        # Data/Pandas stuff
        # ---------------------------

        df = self._to_functional_df(
            annotation=annotation,
            metric=metric,
            taxa_stratified=False,
            fill_missing=True,
        )
        num_of_items = len(df.index)
        ocx_feature_name_map = df.ocx_feature_name_map

        def get_feature_name(feature_id):
            return ocx_feature_name_map.get(feature_id) or feature_id or ""

        # Using fast pandas/numpy approach to select top N functions
        agg_row = df.mean()
        agg_row.sort_values(ascending=False, inplace=True)
        to_keep = agg_row[:num_of_functions]  # Needed to sort y-axis
        to_drop = agg_row[num_of_functions:]

        df.drop(columns=to_drop.index, inplace=True)

        # Joining with metadata
        metadata = df.ocx_metadata
        if metadata.index.name != "sample_id":
            metadata.set_index("sample_id", drop=False, inplace=True)
        metadata.drop("created_at", axis=1, inplace=True)

        # Preparing "Label" column before the merge. So function ids would not override metadata columns
        if label is not None:
            metadata["Label"] = self._make_labels_by_item_id(metadata, label)
        else:
            metadata["Label"] = metadata["sample_name"]
        df = df.join(metadata)

        # Wide-form data -> Long-form data
        df = df.melt(
            id_vars=list(metadata.columns),
            var_name="function_id",
            value_name="value",
        )
        # It is helpful to have function_id and function_name not just one of them
        df["function_name"] = pd.Series([get_feature_name(x) for x in df["function_id"]])

        column_kwargs = {}
        if haxis:
            column_kwargs = {
                "column": alt.Column(
                    haxis,
                    type="nominal",
                    header=alt.Header(titleOrient="bottom", labelOrient="bottom"),
                ),
            }

        # Sorting X/Y axis
        if sort_x:
            sort_x_values = sort_helper(sort_x, df["Label"].tolist())
        else:
            sort_x_values = None

        sort_y_values = [get_feature_name(x) for x in to_keep.index]

        # Altair chart
        # ------------

        chart = (
            alt.Chart(df)
            .mark_rect()
            .encode(
                x=alt.X("Label:N", title=xlabel, sort=sort_x_values),
                y=alt.Y(
                    "function_name:N", title="Function", sort=sort_y_values
                ),
                color=alt.Color("value:Q", title=metric.name),
                tooltip=[
                    alt.Tooltip("Label:N", title="Label"),
                    alt.Tooltip("function_name:N", title="Function Name"),
                    alt.Tooltip("function_id:N", title="Function ID"),
                    alt.Tooltip("value:Q", format=".02f", title=metric.name),
                ],
                **column_kwargs,
            )
        )

        chart = chart.properties(
            **prepare_props(
                title=title,
                width=width or (15 * num_of_items),
                height=height or (15 * num_of_functions),
            )
        )
        if haxis:
            chart = chart.resolve_scale(x="independent")

        if return_chart:
            return chart
        else:
            chart.interactive().display()
