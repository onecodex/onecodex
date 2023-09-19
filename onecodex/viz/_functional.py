import pandas as pd

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

        df.drop(columns=to_drop.index, inplace=True)

        # TODO: comment to explain
        metadata = df.ocx_metadata
        if metadata.index.name != "sample_id":
            metadata.index = pd.Index(metadata["sample_id"], name="sample_id")
        metadata.drop("created_at", axis=1, inplace=True)
        with_metadata = df.join(metadata)

        data = with_metadata.melt(
            id_vars=list(metadata.columns),  # TODO: do we need a functional_id ?
            var_name="function_id",
            value_name="value",
        )
        data["function_id"] = data["function_id"].apply(
            lambda fid: df.ocx_feature_name_map.get(fid, fid)
        )

        # TODO: add sorting. This is default
        y_sort = list(to_keep.index)
        y_sort = [df.ocx_feature_name_map.get(x, x) for x in y_sort]

        chart = (
            alt.Chart(data)
            .mark_rect()
            .encode(
                x=alt.X("name:N", title="Sample Name"),  # TODO: from params
                y=alt.Y("function_id:N", title="Function ID", sort=y_sort),  # TODO: Maybe name?
                color=alt.Color("value:Q", title=metric.name),
                tooltip=[
                    alt.Tooltip("name", title="Sample name"),  # TODO: from params
                    alt.Tooltip("function_id", title="Function ID"),  # TODO: Maybe name?
                    alt.Tooltip("value:Q", format=".02f", title=metric.name),
                ],
            )
        )

        chart = chart.properties(**prepare_props(title=title, width=width, height=height))

        if return_chart:
            return chart
        else:
            chart.interactive().display()
