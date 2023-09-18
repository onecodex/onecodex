from onecodex.lib.enums import (
    FunctionalAnnotations,
    FunctionalAnnotationsMetric,
)
from onecodex.viz._primitives import prepare_props


class VizFunctionalHeatmapMixin(object):
    def plot_functional_heatmap(
        self,
        num_of_functions=10,
        annotation=FunctionalAnnotations.Go,
        metric=None,
        return_chart=False,
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
        return_chart : `bool`, optional
            When True, return an `altair.Chart` object instead of displaying the resulting plot in
            the current notebook.
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

        import altair as alt

        annotation = FunctionalAnnotations(annotation)
        if metric is None:
            if annotation == FunctionalAnnotations.Pathways:
                metric = FunctionalAnnotationsMetric.Abundance
            else:
                metric = FunctionalAnnotationsMetric.Cpm
        metric = FunctionalAnnotationsMetric(metric)

        df = self._to_functional_df(
            annotation=annotation,
            metric=metric,
            taxa_stratified=False,
            fill_missing=True,
        )

        # TODO: comment to explain
        agg_row = df.mean()
        agg_row.sort_values(ascending=False, inplace=True)
        to_keep = agg_row[:num_of_functions]
        to_drop = agg_row[num_of_functions:]

        function_id_to_name = df.ocx_feature_name_map

        df.drop(columns=to_drop.index, inplace=True)

        metadata = df.ocx_metadata
        if metadata.index.name != "sample_id":
            metadata = metadata.reindex(metadata["sample_id"])
        metadata.drop("created_at", axis=1, inplace=True)
        with_metadata = df.join(metadata)

        labels = with_metadata["name"]  # TODO: is that the proper label?
        df.insert(0, "__label", labels)
        # TODO: Maybe just use UUID? df.reset_index(names=["__label"], inplace=True)

        # TODO: add sorting. This is default
        y_sort = list(to_keep.index)

        data = with_metadata.melt(
            id_vars=list(metadata.columns),
            var_name="function_id",
            value_name="value",
        )
        data["function_id"] = data["function_id"].apply(
            lambda fid: function_id_to_name.get(fid, fid)
        )

        chart = (
            alt.Chart(data)
            .mark_rect()
            .encode(
                x=alt.X("name:N", title="TODO"),
                y=alt.Y("function_id:N", title="Function ID", sort=y_sort),
                color=alt.Color("value:Q", title=metric.name),
                tooltip=[
                    alt.Tooltip("name", title="TODO"),
                    alt.Tooltip("function_id", title="Function ID"),
                    alt.Tooltip("value:Q", format=".02f", title=metric.name),
                ],
            )
        )

        chart = chart.properties(**prepare_props(title=title, width=width, height=height))

        if return_chart:
            return chart
        else:
            chart.interactive().display()

    def _to_normalized_functional_df(self, annotation, metric):
        result = self._to_functional_df(
            annotation=annotation,
            metric=metric,
            taxa_stratified=False,
            fill_missing=True,
        )
        result.drop(
            axis=1,
            labels=["UNGROUPED", "UNMAPPED", "UNINTEGRATED"],
            inplace=True,
            errors="ignore",
        )
        return result
