from onecodex.lib.enums import (
    FunctionalAnnotations,
    FunctionalAnnotationsMetric,
)
from onecodex.viz._primitives import (
    prepare_props,
)


class VizFunctionalHeatmapMixin(object):
    def plot_functional_heatmap(
        self,
        return_chart=False,
        num_of_functions=10,
        title=None,
        width=None,
        height=None,
    ):
        """TODO"""
        # TODO: num_of_functions validate???
        # TODO: default_size_kwargs = {"width": "container", "height": "container"}

        # TODO: Placeholder. Remove.

        import altair as alt
        import pandas as pd

        def key_func(key):
            return pd.Index([df[c].mean() for c in key])

        df = self._to_functional_df(
            annotation=FunctionalAnnotations.Go,
            metric=FunctionalAnnotationsMetric.Cpm,
            taxa_stratified=False,
            fill_missing=True,
        )
        df.sort_index(axis=1, key=key_func, ascending=False, inplace=True)
        df.drop(
            axis=1,
            labels=["UNGROUPED", "UNMAPPED", "UNINTEGRATED"],
            inplace=True,
            errors="ignore",
        )  # TODO: Pathways has somethings as well
        df = df.iloc[:, :num_of_functions]

        # TODO: Funtional/Sample labels
        labels = [f"Label {i}" for i in range(len(df))]
        df.insert(1, "__label", labels)

        # TODO: chart/tooltip labels
        chart = (
            alt.Chart(
                df.melt(
                    id_vars=["__label"],
                    var_name="function_id",
                    value_name="value",
                ),
            )
            .mark_rect()
            .encode(
                x=alt.X("__label:N", title="File"),
                y=alt.Y("function_id:N", title="Function ID"),
                color=alt.Color("value:Q"),
                tooltip=[
                    alt.Tooltip("__label", title="File"),
                    alt.Tooltip("function_id", title="Function ID"),
                    alt.Tooltip("value:Q", format=".02f", title="CPM"),
                ],
            )
        )

        chart = chart.properties(**prepare_props(title=title, width=width, height=height))

        if return_chart:
            return chart
        else:
            chart.interactive().display()
