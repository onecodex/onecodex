import warnings
from typing import Optional, Literal, Union

from onecodex.lib.enums import AlphaDiversityMetric, Rank, BaseEnum
from onecodex.exceptions import OneCodexException, PlottingException, PlottingWarning
from onecodex.viz._primitives import (
    prepare_props,
    sort_helper,
    get_classification_url,
    escape_chart_fields,
)


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
        width: Optional[Union[int, Literal["container"]]] = 200,
        height: Optional[Union[int, Literal["container"]]] = 400,
        facet_by=None,
        coerce_haxis_dates=True,
        secondary_haxis=None,
        match_taxonomy=True,
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
            - {'simpson', 'observed_taxa', 'shannon'}: an alpha diversity statistic to calculate for
              each sample. Note that Shannon diversity is calculated using log base ``e`` (natural
              log).
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

        width : `int` or `str`, optional
            Sets `altair.Chart.width`. If `"container"`, chart width will respond to the width of
            the HTML container it is rendered in.

        height : `int` or `str`, optional
            Sets `altair.Chart.height`. If `"container"`, chart height will respond to the height of
            the HTML container it is rendered in.

        facet_by : `string`, optional
            The metadata field used to facet samples by (i.e. to create a separate subplot for each
            group of samples).

        coerce_haxis_dates : `bool`, optional
            If ``True``, ``haxis`` field name(s) containing the word "date" (after splitting on
            underscores) will be coerced to datetime dtype. For example, the field "date_collected"
            will be coerced if ``coerce_haxis_dates`` is ``True``.

        secondary_haxis : str or tuple of str, optional
            The secondary metadata field (or tuple containing multiple categorical fields) to be
            plotted on the horizontal axis.

        match_taxonomy : `bool`, default=True
            Whether or not to consider taxonomic names when looking for metadata fields mapped to
            plot attributes including `vaxis`, `haxis`, `secondary_axis`, `facet_by`, & `label`

        Examples
        --------
        Generate a boxplot of the abundance of Bacteroides (genus) of samples grouped by whether the
        individuals are allergic to dogs, cats, both, or neither.

        >>> plot_metadata(haxis=('allergy_dogs', 'allergy_cats'), vaxis='Bacteroides')
        """
        # Deferred imports
        import altair as alt
        import pandas as pd

        if rank is None:
            raise OneCodexException("Please specify a rank or 'auto' to choose automatically")

        if not PlotType.has_value(plot_type):
            raise OneCodexException("Plot type must be one of: auto, boxplot, scatter")

        if haxis == secondary_haxis:
            raise OneCodexException("`haxis` and `secondary_haxis` cannot be the same field(s).")

        if len(self._results) == 0:
            raise PlottingException(
                "There are too few samples for metadata plots after filtering. Please select 1 or "
                "more samples to plot."
            )

        # alpha diversity is only allowed on vertical axis--horizontal can be magically mapped
        metadata_fields = [haxis, "Label"]
        if facet_by:
            metadata_fields.append(facet_by)
        if secondary_haxis:
            metadata_fields.append(secondary_haxis)

        metadata_results = self._metadata_fetch(
            metadata_fields, label=label, match_taxonomy=match_taxonomy
        )
        df = metadata_results.df
        magic_fields = metadata_results.renamed_fields

        exclude_classifications_without_abundances = True
        if AlphaDiversityMetric.has_value(vaxis):
            df.loc[:, vaxis] = self.alpha_diversity(vaxis, rank=rank)
            magic_fields[vaxis] = vaxis
            df.dropna(subset=[magic_fields[vaxis]], inplace=True)
        else:
            # if it's not alpha diversity, vertical axis can also be magically mapped
            vert_metadata_results = self._metadata_fetch([vaxis], match_taxonomy=match_taxonomy)
            vert_df = vert_metadata_results.df
            vert_magic_fields = vert_metadata_results.renamed_fields
            exclude_classifications_without_abundances = (
                vaxis in vert_metadata_results.taxonomy_fields
            )

            # we require the vertical axis to be numerical otherwise plots get weird
            if (
                pd.api.types.is_bool_dtype(vert_df[vert_magic_fields[vaxis]])
                or isinstance(vert_df[vert_magic_fields[vaxis]].dtype, pd.CategoricalDtype)
                or pd.api.types.is_object_dtype(vert_df[vert_magic_fields[vaxis]])
                or not pd.api.types.is_numeric_dtype(vert_df[vert_magic_fields[vaxis]])
            ):
                raise OneCodexException("Metadata field on vertical axis must be numerical")

            df = pd.concat([df, vert_df], axis=1).dropna(subset=[vert_magic_fields[vaxis]])
            magic_fields.update(vert_magic_fields)

        if (
            exclude_classifications_without_abundances
            and self._classification_ids_without_abundances
        ):
            df = df[~df.index.isin(self._classification_ids_without_abundances)]

            if len(df) == 0:
                raise PlottingException(
                    "Abundances are not calculated for any of the selected samples. Please select "
                    "a different metric or a different set of samples to plot."
                )

            warnings.warn(
                f"{len(self._classification_ids_without_abundances)} sample(s) have no abundances "
                f"calculated and have been omitted from the metadata plot.",
                PlottingWarning,
            )

        # plots can look different depending on what the horizontal axis contains
        if pd.api.types.is_datetime64_any_dtype(df[magic_fields[haxis]]):
            if plot_type == PlotType.Auto:
                plot_type = PlotType.BoxPlot
        elif "date" in magic_fields[haxis].split("_"):
            if coerce_haxis_dates:
                df[magic_fields[haxis]] = pd.to_datetime(df[magic_fields[haxis]], utc=True)
            if plot_type == PlotType.Auto:
                plot_type = PlotType.BoxPlot
        elif (
            pd.api.types.is_bool_dtype(df[magic_fields[haxis]])
            or isinstance(df[magic_fields[haxis]].dtype, pd.CategoricalDtype)
            or pd.api.types.is_object_dtype(df[magic_fields[haxis]])
        ):
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
            if facet_by:
                xlabel = ""
            else:
                xlabel = magic_fields[haxis]

        if ylabel is None:
            ylabel = magic_fields[vaxis]

        encode_kwargs = {}
        if facet_by:
            encode_kwargs["column"] = alt.Column(
                facet_by,
                header=alt.Header(title=facet_by, titleOrient="bottom", labelOrient="bottom"),
            )

        x_kwargs = {"axis": alt.Axis(title=xlabel)}
        if secondary_haxis:
            encode_kwargs["xOffset"] = secondary_haxis
            encode_kwargs["color"] = secondary_haxis

        if plot_type == "scatter":
            df = df.reset_index()
            sort_order = sort_helper(sort_x, df[magic_fields[haxis]].tolist())
            df["url"] = df["classification_id"].apply(get_classification_url)

            encode_kwargs.update(
                dict(
                    x=alt.X(magic_fields[haxis], sort=sort_order, **x_kwargs),
                    y=alt.Y(magic_fields[vaxis], axis=alt.Axis(title=ylabel)),
                    tooltip=["Label", "{}:Q".format(magic_fields[vaxis])],
                    href="url:N",
                )
            )
            chart = alt.Chart(df).mark_circle().encode(**encode_kwargs)
        elif plot_type == PlotType.BoxPlot:
            if sort_x:
                raise OneCodexException("Must not specify sort_x when plot_type is boxplot")

            boxplot_kwargs = {"median": {"stroke": "black"}}
            if not secondary_haxis or width == "container":
                box_size = 45
                increment = 5
                n_boxes = len(df[magic_fields[haxis]].unique())

                if width and width != "container" and (n_boxes * (box_size + increment)) > width:
                    box_size = ((width / n_boxes) // increment) * increment - increment

                boxplot_kwargs["size"] = box_size

            chart = (
                alt.Chart(df)
                .mark_boxplot(**boxplot_kwargs)
                .encode(
                    x=alt.X(magic_fields[haxis], **x_kwargs),
                    y=alt.Y(magic_fields[vaxis], axis=alt.Axis(title=ylabel)),
                    **encode_kwargs,
                )
            )

        if facet_by and not secondary_haxis:
            chart = chart.resolve_scale(x="independent")

        chart = chart.properties(**prepare_props(title=title, height=height, width=width))
        escape_chart_fields(chart)

        if return_chart:
            return chart
        else:
            chart.interactive().display()
