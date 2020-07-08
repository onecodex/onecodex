import warnings

from onecodex.lib.enums import AlphaDiversityMetric, Rank, BaseEnum
from onecodex.exceptions import OneCodexException, PlottingException, PlottingWarning
from onecodex.viz._primitives import prepare_props, sort_helper


class PlotType(BaseEnum):
    Auto = "auto"
    BoxPlot = "boxplot"
    Scatter = "scatter"


class VizMetadataMixin(object):
    def plot_metadata(
        self,
        rank=Rank.Auto,
        haxis="Label",
        vaxis=AlphaDiversityMetric.Shannon,
        title=None,
        xlabel=None,
        ylabel=None,
        return_chart=False,
        plot_type=PlotType.Auto,
        label=None,
        sort_x=None,
        width=200,
        height=400,
    ):
        """Plot an arbitrary metadata field versus an arbitrary quantity as a boxplot or scatter plot.

        Parameters
        ----------
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.

        haxis : `string`, optional
            The metadata field (or tuple containing multiple categorical fields) to be plotted on
            the horizontal axis.

        vaxis : `string`, optional
            Data to be plotted on the vertical axis. Can be any one of the following:

            - A metadata field: the name of a metadata field containing numerical data
            - {'simpson', 'chao1', 'shannon'}: an alpha diversity statistic to calculate for each sample
            - A taxon name: the name of a taxon in the analysis
            - A taxon ID: the ID of a taxon in the analysis

        title : `string`, optional
            Text label at the top of the plot.

        xlabel : `string`, optional
            Text label along the horizontal axis.

        ylabel : `string`, optional
            Text label along the vertical axis.

        plot_type : {'auto', 'boxplot', 'scatter'}
            By default, will determine plot type automatically based on the data. Otherwise, specify
            one of 'boxplot' or 'scatter' to set the type of plot manually.

        label : `string` or `callable`, optional
            A metadata field (or function) used to label each analysis. If passing a function, a
            dict containing the metadata for each analysis is passed as the first and only
            positional argument. The callable function must return a string.

        sort_x : `list` or `callable`, optional
            Either a list of sorted labels or a function that will be called with a list of x-axis labels
            as the only argument, and must return the same list in a user-specified order.

        Examples
        --------
        Generate a boxplot of the abundance of Bacteroides (genus) of samples grouped by whether the
        individuals are allergy to dogs, cats, both, or neither.

        >>> plot_metadata(haxis=('allergy_dogs', 'allergy_cats'), vaxis='Bacteroides')
        """
        # Deferred imports
        import altair as alt
        import pandas as pd

        if rank is None:
            raise OneCodexException("Please specify a rank or 'auto' to choose automatically")

        if not PlotType.has_value(plot_type):
            raise OneCodexException("Plot type must be one of: auto, boxplot, scatter")

        if len(self._results) < 1:
            raise PlottingException(
                "There are too few samples for metadata plots after filtering. Please select 1 or "
                "more samples to plot."
            )

        # alpha diversity is only allowed on vertical axis--horizontal can be magically mapped
        df, magic_fields = self._metadata_fetch([haxis, "Label"], label=label)

        if AlphaDiversityMetric.has_value(vaxis):
            df.loc[:, vaxis] = self.alpha_diversity(vaxis, rank=rank)
            magic_fields[vaxis] = vaxis
        else:
            # if it's not alpha diversity, vertical axis can also be magically mapped
            vert_df, vert_magic_fields = self._metadata_fetch([vaxis])

            # we require the vertical axis to be numerical otherwise plots get weird
            if (
                pd.api.types.is_bool_dtype(vert_df[vert_magic_fields[vaxis]])
                or pd.api.types.is_categorical_dtype(vert_df[vert_magic_fields[vaxis]])
                or pd.api.types.is_object_dtype(vert_df[vert_magic_fields[vaxis]])
                or not pd.api.types.is_numeric_dtype(vert_df[vert_magic_fields[vaxis]])
            ):  # noqa
                raise OneCodexException("Metadata field on vertical axis must be numerical")

            df = pd.concat([df, vert_df], axis=1).dropna(subset=[vert_magic_fields[vaxis]])
            magic_fields.update(vert_magic_fields)

        # plots can look different depending on what the horizontal axis contains
        if pd.api.types.is_datetime64_any_dtype(df[magic_fields[haxis]]):
            if plot_type == PlotType.Auto:
                plot_type = PlotType.BoxPlot
        elif "date" in magic_fields[haxis].split("_"):
            df.loc[:, magic_fields[haxis]] = df.loc[:, magic_fields[haxis]].apply(
                pd.to_datetime, utc=True
            )

            if plot_type == PlotType.Auto:
                plot_type = PlotType.BoxPlot
        elif (
            pd.api.types.is_bool_dtype(df[magic_fields[haxis]])
            or pd.api.types.is_categorical_dtype(df[magic_fields[haxis]])
            or pd.api.types.is_object_dtype(df[magic_fields[haxis]])
        ):  # noqa
            df = df.fillna({field: "N/A" for field in df.columns})

            if plot_type == PlotType.Auto:
                # if data is categorical but there is only one value per sample, scatter plot instead
                if len(df[magic_fields[haxis]].unique()) == len(df[magic_fields[haxis]]):
                    plot_type = PlotType.Scatter
                else:
                    plot_type = PlotType.BoxPlot
        elif pd.api.types.is_numeric_dtype(df[magic_fields[haxis]]):
            df = df.dropna(subset=[magic_fields[vaxis]])

            if plot_type == PlotType.Auto:
                plot_type = PlotType.Scatter
        else:
            raise OneCodexException(
                "Unplottable column type for horizontal axis ({})".format(haxis)
            )

        if xlabel is None:
            xlabel = magic_fields[haxis]

        if ylabel is None:
            ylabel = magic_fields[vaxis]

        if plot_type == "scatter":
            df = df.reset_index()

            sort_order = sort_helper(sort_x, df[magic_fields[haxis]].tolist())

            alt_kwargs = dict(
                x=alt.X(magic_fields[haxis], axis=alt.Axis(title=xlabel), sort=sort_order),
                y=alt.Y(magic_fields[vaxis], axis=alt.Axis(title=ylabel)),
                tooltip=["Label", "{}:Q".format(magic_fields[vaxis])],
                href="url:N",
                url="https://app.onecodex.com/classification/" + alt.datum.classification_id,
            )

            chart = (
                alt.Chart(df)
                .transform_calculate(url=alt_kwargs.pop("url"))
                .mark_circle()
                .encode(**alt_kwargs)
            )

        elif plot_type == PlotType.BoxPlot:
            if sort_x:
                raise OneCodexException("Must not specify sort_x when plot_type is boxplot")

            # See the following issue in case this gets fixed in altair:
            # https://github.com/altair-viz/altair/issues/2144
            if (df.groupby(magic_fields[haxis]).size() < 2).any():
                warnings.warn(
                    "There is at least one sample group consisting of only a single sample. Groups "
                    "of size 1 may not have their boxes displayed in the plot.",
                    PlottingWarning,
                )

            chart = (
                alt.Chart(df)
                .mark_boxplot(size=35)
                .encode(
                    x=alt.X(magic_fields[haxis], axis=alt.Axis(title=xlabel)),
                    y=alt.Y(magic_fields[vaxis], axis=alt.Axis(title=ylabel)),
                )
            )

        chart = chart.properties(**prepare_props(title=title, height=height, width=width))

        if return_chart:
            return chart
        else:
            chart.interactive().display()
