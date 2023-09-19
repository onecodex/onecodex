from onecodex.lib.enums import (
    FunctionalAnnotations,
    FunctionalAnnotationsMetric,
)
from onecodex.viz._primitives import prepare_props, sort_helper, get_analysis_url


class VizFunctionalHeatmapMixin(object):
    def plot_functional_heatmap(
        self,
        num_of_functions=10,
        annotation=FunctionalAnnotations.Go,
        metric=None,
        sort_x=None,
        label=None,
        return_chart=False,
        xlabel="",
        title=None,
        width=None,
        height=None,
    ):
        """Plot a TODO heatmap.

        Parameters
        ----------
        num_of_functions : `int`, optional
            TODO
        annotation : `FunctionalAnnotations` or `str`, optional
            TODO
        metric : `FunctionalAnnotationsMetric` or `str`, optional
            TODO
        sort_x : `list` or `callable`, optional
            Either a list of sorted labels or a function that will be called with a list of x-axis labels
            as the only argument, and must return the same list in a user-specified order.
        label : `str` or `callable`, optional
            A metadata field (or function) used to label each analysis. If passing a function, a
            dict containing the metadata for each analysis is passed as the first and only
            positional argument. The callable function must return a string.
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

        Examples
        --------
        TODO

        >>> "TODO"
        """
        # TODO: num_of_functions validate???

        # Deferred imports
        import altair as alt

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
        ocx_feature_name_map = df.ocx_feature_name_map

        # TODO: comment to explain
        agg_row = df.drop("functional_profile_id", axis=1).mean()
        agg_row.sort_values(ascending=False, inplace=True)
        to_keep = agg_row[:num_of_functions]
        to_drop = agg_row[num_of_functions:]

        df.drop(columns=to_drop.index, inplace=True)

        # TODO: comment to explain
        metadata = df.ocx_metadata
        if metadata.index.name != "sample_id":
            metadata.set_index("sample_id", drop=False, inplace=True)
        metadata.drop("created_at", axis=1, inplace=True)
        if label is not None:
            # TODO: double check if other plots are ignoring `Labels` and just overriding it.
            metadata["Label"] = self._make_labels_by_item_id(metadata, label)
        else:
            metadata["Label"] = metadata["name"]  # TODO: can we assume that name exists?

        df["url"] = df["functional_profile_id"].apply(lambda fpid: get_analysis_url(fpid))
        df = df.join(metadata)

        # Wide-form data -> Long-form data
        df = df.melt(
            id_vars=list(metadata.columns) + ["functional_profile_id", "url"],
            var_name="function_id",
            value_name="value",
        )
        df["function_id"] = df["function_id"].apply(lambda fid: ocx_feature_name_map.get(fid, fid))

        # Sorting X/Y axis
        if sort_x:
            sort_x_values = sort_helper(sort_x, df["Label"].tolist())
        else:
            sort_x_values = None

        sort_y_values = list(to_keep.index)
        sort_y_values = [ocx_feature_name_map.get(x, x) for x in sort_y_values]

        # Altair chart
        # ------------

        chart = (
            alt.Chart(df)
            .mark_rect()
            .encode(
                x=alt.X("Label:N", title=xlabel, sort=sort_x_values),
                y=alt.Y(
                    "function_id:N", title="Function ID", sort=sort_y_values
                ),  # TODO: Maybe name?
                color=alt.Color("value:Q", title=metric.name),
                tooltip=[
                    alt.Tooltip("Label", title="Label"),  # TODO: maybe change title ?
                    # TODO: Display function ID **and** name?
                    alt.Tooltip("function_id", title="Function ID"),
                    alt.Tooltip("value:Q", format=".02f", title=metric.name),
                ],
                href="url:N",
            )
        )

        chart = chart.properties(**prepare_props(title=title, width=width, height=height))

        if return_chart:
            return chart
        else:
            chart.interactive().display()
