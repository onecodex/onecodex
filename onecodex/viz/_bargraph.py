import altair as alt
from onecodex.exceptions import OneCodexException


class VizBargraphMixin(object):
    def plot_bargraph(
        self,
        rank="auto",
        normalize="auto",
        top_n="auto",
        threshold="auto",
        title=None,
        xlabel=None,
        ylabel=None,
        tooltip=None,
        return_chart=False,
        legend="auto",
    ):
        """Plot a bargraph of relative abundance of taxa for multiple samples.

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

        Examples
        --------
        Plot a bargraph of the top 10 most abundant genera

        >>> plot_bargraph(rank='genus', top_n=10)
        """
        if rank is None:
            raise OneCodexException(
                "Please specify a rank or 'auto' to choose automatically"
            )

        if not (threshold or top_n):
            raise OneCodexException("Please specify at least one of: threshold, top_n")

        if top_n == "auto" and threshold == "auto":
            top_n = 10
            threshold = None
        elif top_n == "auto" and threshold != "auto":
            top_n = None
        elif top_n != "auto" and threshold == "auto":
            threshold = None

        if legend == 'auto':
            legend = self._field

        df = self.to_df(
            rank=rank,
            normalize=normalize,
            top_n=top_n,
            threshold=threshold,
            table_format="long",
        )

        if tooltip:
            if not isinstance(tooltip, list):
                tooltip = [tooltip]
        else:
            tooltip = []

        tooltip.append("Label")

        # takes metadata columns and returns a dataframe with just those columns
        # renames columns in the case where columns are taxids
        magic_metadata, magic_fields = self._metadata_fetch(tooltip)

        for f in tooltip:
            df[magic_fields[f]] = magic_metadata[magic_fields[f]][
                df["classification_id"]
            ].tolist()

        # add taxa names
        df["tax_name"] = [
            "{} ({})".format(self.taxonomy["name"][t], t)
            if t in self.taxonomy["name"]
            else t
            for t in df["tax_id"]
        ]

        #
        # TODO: how to sort bars in bargraph
        # - abundance (mean across all samples)
        # - parent taxon (this will require that we make a few assumptions
        # about taxonomic ranks but as all taxonomic data will be coming from
        # OCX this should be okay)
        #

        ylabel = self._field if ylabel is None else ylabel
        xlabel = 'Label' if xlabel is None else xlabel

        chart = (
            alt.Chart(df)
            .mark_bar()
            .encode(
                x=alt.X("Label", axis=alt.Axis(title=xlabel)),
                y=alt.Y(self._field, axis=alt.Axis(title=ylabel)),
                color=alt.Color("tax_name", legend=alt.Legend(title=legend)),
                tooltip=["{}:Q".format(self._field)] + [magic_fields[f] for f in tooltip],
                href="url:N",
            )
        )

        # add title to chart
        # (cannot specify None or False for no title)
        chart = chart.properties(title=title) if title else chart

        return chart if return_chart else chart.display()
