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
        haxis=None,
        legend="auto",
        label=None,
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
        haxis : `string`, optional
            The metadata field (or tuple containing multiple categorical fields) used to group
            samples together.
        legend: `string`, optional
            Title for color scale. Defaults to the field used to generate the plot, e.g.
            readcount_w_children or abundance.
        label : `string` or `callable`, optional
            A metadata field (or function) used to label each analysis. If passing a function, a
            dict containing the metadata for each analysis is passed as the first and only
            positional argument. The callable function must return a string.

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

        if haxis:
            tooltip.append(haxis)

        tooltip.insert(0, "Label")

        # takes metadata columns and returns a dataframe with just those columns
        # renames columns in the case where columns are taxids
        magic_metadata, magic_fields = self._metadata_fetch(tooltip, label=label)

        # add sort order to long-format df
        if haxis:
            sort_order = magic_metadata.sort_values(magic_fields[haxis]).index.tolist()

            for sort_num, sort_class_id in enumerate(sort_order):
                magic_metadata.loc[sort_class_id, 'sort_order'] = sort_num

            df['sort_order'] = magic_metadata['sort_order'][
                df["classification_id"]
            ].tolist()

            sort_order = alt.EncodingSortField(field='sort_order', op="mean")
        else:
            sort_order = None

        # transfer metadata from wide-format df (magic_metadata) to long-format df
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
        xlabel = '' if xlabel is None else xlabel

        # should ultimately be Label, tax_name, readcount_w_children, then custom fields
        tooltip_for_altair = [magic_fields[f] for f in tooltip]
        tooltip_for_altair.insert(1, "tax_name")
        tooltip_for_altair.insert(2, "{}:Q".format(self._field))

        # generate dataframes to plot, one per facet
        dfs_to_plot = []

        if haxis:
            # if using facets, first facet is just the vertical axis
            blank_df = df.iloc[:1].copy()
            blank_df[self._field] = 0

            dfs_to_plot.append(blank_df)

            for md_val in magic_metadata[magic_fields[haxis]].unique():
                plot_df = df.where(df[magic_fields[haxis]] == md_val).dropna()

                # preserve booleans
                if magic_metadata[magic_fields[haxis]].dtype == 'bool':
                    plot_df[magic_fields[haxis]] = plot_df[magic_fields[haxis]].astype(bool)

                dfs_to_plot.append(plot_df)
        else:
            dfs_to_plot.append(df)

        charts = []

        for plot_num, plot_df in enumerate(dfs_to_plot):
            chart = (
                alt.Chart(plot_df)
                .mark_bar()
                .encode(
                    x=alt.X("Label", axis=alt.Axis(title=xlabel), sort=sort_order),
                    y=alt.Y(self._field, axis=alt.Axis(title=ylabel), scale=alt.Scale(domain=[0, 1], zero=True, nice=False)),
                    color=alt.Color("tax_name", legend=alt.Legend(title=legend)),
                    tooltip=tooltip_for_altair,
                    href="url:N",
                )
            )

            if haxis:
                if plot_num == 0:
                    # first plot (blank_df) has vert axis but no horiz axis
                    chart.encoding.x.axis = None
                elif plot_num > 0:
                    # strip vertical axis from subsequent facets
                    chart.encoding.y.axis = None

                    # facet's title set to value of metadata in this group
                    chart.title = str(plot_df[magic_fields[haxis]].tolist()[0])

            charts.append(chart)

        # add all the facets together
        final_chart = charts[0]

        if len(charts) > 1:
            for chart in charts[1:]:
                final_chart |= chart

        # add title to chart
        # (cannot specify None or False for no title)
        final_chart = final_chart.properties(title=title) if title else final_chart
        return final_chart if return_chart else final_chart.display()
