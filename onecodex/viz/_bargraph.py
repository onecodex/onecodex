from onecodex.exceptions import OneCodexException, PlottingException
from onecodex.lib.enums import Link, Metric, Rank
from onecodex.viz._primitives import (
    get_classification_url,
    get_ncbi_taxonomy_browser_url,
    get_unique_column,
    interleave_palette,
    prepare_props,
    sort_helper,
    escape_chart_fields,
)


class VizBargraphMixin(object):
    def plot_bargraph(
        self,
        rank=Rank.Auto,
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
        sort_x=None,
        include_taxa_missing_rank=None,
        include_other=True,
        width=None,
        height=None,
        group_by=None,
        link=Link.Ocx,
        match_taxonomy=True,
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
            The metadata field (or tuple containing multiple categorical fields) used to facet
            samples.
        legend: `string` or `altair.Legend`, optional
            If a string is provided, it will be used as the legend title. Defaults to the metric
            used to generate the plot, e.g. readcount_w_children or abundance. Alternatively, an
            `altair.Legend` instance may be provided for legend customization.
        label : `string` or `callable`, optional
            A metadata field (or function) used to label each analysis. If passing a function, a
            dict containing the metadata for each analysis is passed as the first and only
            positional argument. The callable function must return a string.
        sort_x : `list` or `callable`, optional
            Either a list of sorted labels or a function that will be called with a list of x-axis labels
            as the only argument, and must return the same list in a user-specified order.
        include_taxa_missing_rank : `bool`, optional
            Whether or not a row should be plotted for taxa that do not have a designated parent at `rank`.
        group_by : `string`, optional
            The metadata field used to group samples together. Readcounts or abundances will be
            averaged within each group.
        link: {'ocx', 'ncbi'}, optional
            If `link` is 'ocx', clicking a sample will open its classification results in the One
            Codex app. If `link` is 'ncbi', clicking a taxon will open the NCBI taxonomy browser.
        match_taxonomy : `bool`, default=True
            Whether or not to consider taxonomic names when looking for metadata fields mapped to
            plot attributes including `group_by`, `label`, `haxis`, & `tooltip`

        Examples
        --------
        Plot a bargraph of the top 10 most abundant genera

        >>> samples.plot_bargraph(rank='genus', top_n=10)
        """

        # Deferred imports
        import altair as alt

        if rank is None:
            raise OneCodexException("Please specify a rank or 'auto' to choose automatically")

        if not (threshold or top_n):
            raise OneCodexException("Please specify at least one of: threshold, top_n")

        if len(self._results) < 1:
            raise PlottingException(
                "There are too few samples for bargraph plots after filtering. Please select 1 or "
                "more samples to plot."
            )

        if not normalize and self._guess_normalized():
            raise OneCodexException("Data has already been normalized and this cannot be undone.")

        if group_by:
            if not all(kwarg is None for kwarg in (tooltip, haxis, label, sort_x)):
                raise OneCodexException(
                    "`tooltip`, `haxis`, `label`, and `sort_x` are not supported with `group_by`."
                )
            if group_by not in self.metadata:
                raise OneCodexException(
                    f"Metadata field {group_by} not found. Choose from: {', '.join(self.metadata.keys())}"
                )
            if (
                self._metric in {Metric.Readcount, Metric.ReadcountWChildren}
                and self._guess_normalized()
            ):
                raise OneCodexException(
                    "`group_by` may not be used with readcounts that have already been normalized."
                )

        if include_taxa_missing_rank is None:
            include_taxa_missing_rank = self._metric in {
                Metric.ReadcountWChildren,
                Metric.AbundanceWChildren,
            }

        # We're intentionally *not* normalizing or filtering by top_n/threshold in to_df() in case
        # group_by was passed. Grouping samples needs to happen *before* normalization.
        df = self.to_df(
            rank=rank,
            top_n=None,
            threshold=None,
            normalize=None,
            include_taxa_missing_rank=include_taxa_missing_rank,
        )
        pretty_metric_name = self.metric

        if group_by:
            df = df.fillna(0.0).join(self.metadata[group_by]).groupby(group_by, dropna=False).mean()
            # Nicer display for missing metadata values than `null`
            df.index = df.index.fillna("N/A")
            pretty_metric_name = f"Mean {pretty_metric_name}"

        if normalize and (not self._guess_normalized() or group_by):
            # Replace nans with zeros for samples that have a total abundance of zero.
            df = df.div(df.sum(axis=1), axis=0).fillna(0.0)

        # Keep track of empty rows *before* filtering taxa by threshold/top_n.
        # We'll use this below to calculate "Other".
        empty_rows = df[df.sum(axis=1) == 0.0].index

        if top_n == "auto" and threshold == "auto":
            top_n = 10
            threshold = None
        elif top_n == "auto" and threshold != "auto":
            top_n = None
        elif top_n != "auto" and threshold == "auto":
            threshold = None

        if threshold:
            df = df.loc[:, df.max() >= threshold]

        if top_n:
            df = df.loc[:, df.mean().sort_values(ascending=False).iloc[:top_n].index]

        if include_other and normalize:
            # if there are no abundances in the dataframe, df.apply will yield
            # a single item that has `None` as its name. Therefore, we must
            # check if row.name is None AND whether or not the row name is
            # in empty_rows
            df["Other"] = df.apply(
                lambda row: 0.0 if (row.name is None or row.name in empty_rows) else 1 - row.sum(),
                axis=1,
            )

        if isinstance(legend, str):
            if legend == "auto":
                legend = pretty_metric_name
            legend = alt.Legend(title=legend, symbolLimit=40, labelLimit=0)

        if not isinstance(legend, alt.Legend):
            raise TypeError(f"`legend` must be of type str or altair.Legend, not {type(legend)}")

        if tooltip:
            if isinstance(tooltip, list):
                tooltip = tooltip.copy()
            else:
                tooltip = [tooltip]
        else:
            tooltip = []

        if group_by:
            tooltip.append(group_by)
            metadata_columns = []
        else:
            tooltip.insert(0, "Label")
            if haxis:
                tooltip.append(haxis)

            # takes metadata columns and returns a dataframe with just those columns
            # renames columns in the case where columns are taxids
            metadata_results = self._metadata_fetch(
                tooltip, label=label, match_taxonomy=match_taxonomy
            )
            magic_metadata = metadata_results.df
            magic_fields = metadata_results.renamed_fields
            df = df.join(magic_metadata)
            metadata_columns = magic_metadata.columns.tolist()
            tooltip = [magic_fields[f] for f in tooltip]

        id_vars = [df.index.name] + metadata_columns
        tax_id_column = get_unique_column("tax_id", id_vars)
        tax_name_column = get_unique_column("tax_name", id_vars)

        # should ultimately be Label/`group_by`, tax_name, metric name, then custom fields
        tooltip.insert(1, tax_name_column)
        tooltip.insert(2, "{}:Q".format(pretty_metric_name))

        df = df.reset_index().melt(
            id_vars=id_vars,
            var_name=tax_id_column,
            value_name=pretty_metric_name,
        )

        # add taxa names
        df[tax_name_column] = df[tax_id_column].apply(
            lambda t: "{} ({})".format(self.taxonomy["name"][t], t)
            if t in self.taxonomy["name"]
            else t
        )

        #
        # TODO: how to sort bars in bargraph
        # - abundance (mean across all samples)
        # - parent taxon (this will require that we make a few assumptions
        # about taxonomic ranks but as all taxonomic data will be coming from
        # OCX this should be okay)
        #

        ylabel = ylabel or pretty_metric_name

        if xlabel is None:
            xlabel = group_by if group_by else ""

        encode_kwargs = {}
        if haxis:
            encode_kwargs["column"] = alt.Column(
                haxis, header=alt.Header(title=haxis, titleOrient="bottom", labelOrient="bottom")
            )
        domain = sorted(df[tax_name_column].unique())
        no_level_name = "No {}".format(rank)
        color_range = interleave_palette(set(domain) - {"Other", no_level_name})
        other_color = ["#DCE0E5"]
        no_level_color = ["#eeefe1"]

        if include_taxa_missing_rank and no_level_name in domain:
            domain.remove(no_level_name)
            domain = [no_level_name] + domain
            color_range = no_level_color + color_range

        if include_other and "Other" in domain:
            domain.remove("Other")
            domain = ["Other"] + domain
            color_range = other_color + color_range

        sort_order = None
        if not group_by:
            sort_order = sort_helper(sort_x, df["Label"].tolist())

        df["order"] = df[tax_name_column].apply(domain.index)

        if link == Link.Ocx and not group_by:
            df["url"] = df["classification_id"].apply(get_classification_url)
            encode_kwargs["href"] = "url:N"
        elif link == Link.Ncbi:
            df["url"] = df[tax_id_column].apply(get_ncbi_taxonomy_browser_url)
            encode_kwargs["href"] = "url:N"

        y_scale_kwargs = {"zero": True, "nice": False}
        if normalize:
            y_scale_kwargs["domain"] = [0, 1]

        chart = (
            alt.Chart(df)
            .mark_bar()
            .encode(
                x=alt.X(
                    group_by if group_by else "Label", axis=alt.Axis(title=xlabel), sort=sort_order
                ),
                y=alt.Y(
                    pretty_metric_name,
                    axis=alt.Axis(title=ylabel),
                    scale=alt.Scale(**y_scale_kwargs),
                ),
                color=alt.Color(
                    tax_name_column,
                    legend=legend,
                    sort=domain,
                    scale=alt.Scale(domain=domain, range=color_range),
                ),
                tooltip=tooltip,
                order=alt.Order("order", sort="descending"),
                **encode_kwargs,
            )
        )

        if haxis:
            chart = chart.resolve_scale(x="independent")

        chart = chart.properties(**prepare_props(title=title, width=width, height=height))
        escape_chart_fields(chart)
        return chart if return_chart else chart.display()
