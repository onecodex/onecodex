import altair as alt
import pandas as pd

from onecodex.exceptions import OneCodexException
from onecodex.viz import boxplot


class VizMetadataMixin(object):
    def plot_metadata(
        self,
        rank="auto",
        haxis="Label",
        vaxis="simpson",
        title=None,
        xlabel=None,
        ylabel=None,
        return_chart=False,
        plot_type="auto",
        label=None,
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

        Examples
        --------
        Generate a boxplot of the abundance of Bacteroides (genus) of samples grouped by whether the
        individuals are allergy to dogs, cats, both, or neither.

        >>> plot_metadata(haxis=('allergy_dogs', 'allergy_cats'), vaxis='Bacteroides')
        """
        if rank is None:
            raise OneCodexException("Please specify a rank or 'auto' to choose automatically")

        if plot_type not in ("auto", "boxplot", "scatter"):
            raise OneCodexException("Plot type must be one of: auto, boxplot, scatter")

        # alpha diversity is only allowed on vertical axis--horizontal can be magically mapped
        df, magic_fields = self._metadata_fetch([haxis, "Label"], label=label)

        if vaxis in ("simpson", "chao1", "shannon"):
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
            category_type = "T"

            if plot_type == "auto":
                plot_type = "boxplot"
        elif "date" in magic_fields[haxis].split("_"):
            df.loc[:, magic_fields[haxis]] = df.loc[:, magic_fields[haxis]].apply(
                pd.to_datetime, utc=True
            )

            category_type = "T"

            if plot_type == "auto":
                plot_type = "boxplot"
        elif (
            pd.api.types.is_bool_dtype(df[magic_fields[haxis]])
            or pd.api.types.is_categorical_dtype(df[magic_fields[haxis]])
            or pd.api.types.is_object_dtype(df[magic_fields[haxis]])
        ):  # noqa
            df = df.fillna({field: "N/A" for field in df.columns})

            category_type = "N"

            if plot_type == "auto":
                # if data is categorical but there is only one value per sample, scatter plot instead
                if len(df[magic_fields[haxis]].unique()) == len(df[magic_fields[haxis]]):
                    plot_type = "scatter"
                else:
                    plot_type = "boxplot"
        elif pd.api.types.is_numeric_dtype(df[magic_fields[haxis]]):
            df = df.dropna(subset=[magic_fields[vaxis]])

            category_type = "O"

            if plot_type == "auto":
                plot_type = "scatter"
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

            alt_kwargs = dict(
                x=alt.X(magic_fields[haxis], axis=alt.Axis(title=xlabel)),
                y=alt.Y(magic_fields[vaxis], axis=alt.Axis(title=ylabel)),
                tooltip=["Label", "{}:Q".format(vaxis)],
                href="url:N",
                url="https://app.onecodex.com/classification/" + alt.datum.classification_id,
            )

            chart = (
                alt.Chart(df)
                .transform_calculate(url=alt_kwargs.pop("url"))
                .mark_circle()
                .encode(**alt_kwargs)
            )

            if title:
                chart = chart.properties(title=title)
        elif plot_type == "boxplot":
            chart = boxplot(
                df,
                magic_fields[haxis],
                magic_fields[vaxis],
                category_type=category_type,
                title=title,
                xlabel=xlabel,
                ylabel=ylabel,
            )

        if return_chart:
            return chart
        else:
            chart.interactive().display()
