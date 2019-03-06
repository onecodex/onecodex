import altair as alt
import pandas as pd

from onecodex.exceptions import OneCodexException


class VizHeatmapMixin(object):
    def plot_heatmap(self, rank='auto', normalize='auto', top_n='auto', threshold='auto',
                     title=None, xlabel=None, ylabel=None, tooltip=None, return_chart=False,
                     linkage='average', haxis=None, metric='euclidean', legend='auto',
                     label=None):
        """Plot heatmap of taxa abundance/count data for several samples.

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
        haxis : `string`, optional
            The metadata field (or tuple containing multiple categorical fields) used to group
            samples together. Each group of samples will be clustered independently.
        metric : {'braycurtis', 'manhattan', 'jaccard', 'unifrac', 'unweighted_unifrac}, optional
            Function to use when calculating the distance between two samples.
        linkage : {'average', 'single', 'complete', 'weighted', 'centroid', 'median'}
            The type of linkage to use when clustering axes.
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
        label : `string` or `callable`, optional
            A metadata field (or function) used to label each analysis. If passing a function, a
            dict containing the metadata for each analysis is passed as the first and only
            positional argument. The callable function must return a string.

        Examples
        --------
        Plot a heatmap of the relative abundances of the top 10 most abundant families.

        >>> plot_heatmap(rank='family', top_n=10)
        """
        if rank is None:
            raise OneCodexException('Please specify a rank or \'auto\' to choose automatically')

        if not (threshold or top_n):
            raise OneCodexException('Please specify at least one of: threshold, top_n')

        if len(self._results) < 2:
            raise OneCodexException('`plot_heatmap` requires 2 or more valid classification results.')

        if top_n == 'auto' and threshold == 'auto':
            top_n = 10
            threshold = None
        elif top_n == 'auto' and threshold != 'auto':
            top_n = None
        elif top_n != 'auto' and threshold == 'auto':
            threshold = None

        if legend == 'auto':
            legend = self._field

        df = self.to_df(
            rank=rank,
            normalize=normalize,
            top_n=top_n,
            threshold=threshold,
            table_format='long'
        )

        if tooltip:
            if not isinstance(tooltip, list):
                tooltip = [tooltip]
        else:
            tooltip = []

        if haxis:
            tooltip.append(haxis)

        tooltip.insert(0, "Label")

        magic_metadata, magic_fields = self._metadata_fetch(tooltip, label=label)

        # add columns for prettier display
        df['Label'] = magic_metadata['Label'][df['classification_id']].tolist()
        df['tax_name'] = ['{} ({})'.format(self.taxonomy['name'][t], t) for t in df['tax_id']]

        # and for metadata
        for f in tooltip:
            df[magic_fields[f]] = magic_metadata[magic_fields[f]][df['classification_id']].tolist()

        # if we've already been normalized, we must cluster samples by euclidean distance. beta
        # diversity measures won't work with normalized distances.
        if self._guess_normalized():
            if metric != 'euclidean':
                raise OneCodexException('Results are normalized. Please re-run with metric=euclidean')

            df_sample_cluster = self.to_df(
                rank=rank, normalize=normalize, top_n=top_n, threshold=threshold
            )
            df_taxa_cluster = df_sample_cluster
        else:
            df_sample_cluster = self.to_df(
                rank=rank, normalize=False, top_n=top_n, threshold=threshold
            )

            df_taxa_cluster = self.to_df(
                rank=rank, normalize=normalize, top_n=top_n, threshold=threshold
            )

        if haxis is None:
            # cluster only once
            sample_cluster = df_sample_cluster.ocx._cluster_by_sample(rank=rank, metric=metric, linkage=linkage)
            taxa_cluster = df_taxa_cluster.ocx._cluster_by_taxa(linkage=linkage)

            labels_in_order = magic_metadata['Label'][sample_cluster['ids_in_order']].tolist()
        else:
            if not (pd.api.types.is_bool_dtype(df[magic_fields[haxis]]) or  # noqa
                    pd.api.types.is_categorical_dtype(df[magic_fields[haxis]]) or  # noqa
                    pd.api.types.is_object_dtype(df[magic_fields[haxis]])):  # noqa
                raise OneCodexException('Metadata field on horizontal axis can not be numerical')

            # taxa clustered only once
            taxa_cluster = df_taxa_cluster.ocx._cluster_by_taxa(linkage=linkage)

            # cluster samples for every group of metadata
            groups = magic_metadata[magic_fields[haxis]].unique()
            cluster_by_group = {}

            labels_in_order = []

            plot_data = {'x': [], 'y': [], 'o': [], 'b': []}
            label_data = {'x': [], 'y': [], 'label': []}

            for idx, group in enumerate(groups):
                # if value of metadata field is 'null', we have to use pd.isnull, can't use 'is None'
                if pd.isnull(group):
                    c_ids_in_group = magic_metadata.index[pd.isnull(magic_metadata[magic_fields[haxis]])]
                else:
                    c_ids_in_group = magic_metadata.index[magic_metadata[magic_fields[haxis]] == group]

                if len(c_ids_in_group) == 0:
                    continue

                sample_slice = df_sample_cluster.loc[c_ids_in_group]

                if len(c_ids_in_group) < 3:
                    # clustering not possible in this case
                    cluster_by_group[group] = {
                        'ids_in_order': c_ids_in_group
                    }
                else:
                    cluster_by_group[group] = sample_slice.ocx._cluster_by_sample(
                        rank=rank, metric=metric, linkage=linkage
                    )

                plot_data['x'].append(len(labels_in_order) + 0.25)

                labels_in_order.extend(magic_metadata['Label'][cluster_by_group[group]['ids_in_order']].tolist())

                plot_data['x'].append(len(labels_in_order) - 0.25)
                plot_data['y'].extend([0, 0])
                plot_data['o'].extend([0, 1])
                plot_data['b'].extend([idx, idx])

                label_data['x'].append(sum(plot_data['x'][-2:]) / 2)
                label_data['y'].append(1)
                label_data['label'].append(str(group))

            label_bars = alt.Chart(
                pd.DataFrame(plot_data), width=15 * len(df_sample_cluster.index), height=10
            ).mark_line(
                point=False, opacity=0.5
            ).encode(
                x=alt.X(
                    'x',
                    axis=None,
                    scale=alt.Scale(
                        domain=[0, len(df_sample_cluster.index)],
                        zero=True,
                        nice=False
                    )
                ),
                y=alt.Y(
                    'y',
                    axis=None
                ),
                order='o',
                color=alt.Color(
                    'b:N',
                    scale=alt.Scale(
                        domain=list(range(idx + 1)),
                        range=['black'] * (idx + 1)
                    ),
                    legend=None
                )
            )

            label_text = alt.Chart(
                pd.DataFrame(label_data), width=15 * len(df_sample_cluster.index), height=10
            ).mark_text(
                align='center', baseline='middle'
            ).encode(
                x=alt.X(
                    'x',
                    axis=None,
                    scale=alt.Scale(
                        domain=[0, len(df_sample_cluster.index)],
                        zero=True,
                        nice=False
                    )
                ),
                y=alt.Y(
                    'y',
                    axis=alt.Axis(
                        title=haxis, ticks=False, domain=False, labels=False
                    )),
                text='label'
            )

            top_label = alt.layer(label_text, label_bars)

        # should ultimately be Label, tax_name, readcount_w_children, then custom fields
        tooltip_for_altair = [magic_fields[f] for f in tooltip]
        tooltip_for_altair.insert(1, "tax_name")
        tooltip_for_altair.insert(2, "{}:Q".format(self._field))

        alt_kwargs = dict(
            x=alt.X('Label:N', axis=alt.Axis(title=xlabel), sort=labels_in_order),
            y=alt.Y('tax_name:N', axis=alt.Axis(title=ylabel), sort=taxa_cluster['labels_in_order']),
            color=alt.Color('{}:Q'.format(self._field), legend=alt.Legend(title=legend)),
            tooltip=tooltip_for_altair,
            href='url:N',
            url='https://app.onecodex.com/classification/' + alt.datum.classification_id
        )

        chart = alt.Chart(
            df,
            width=15 * len(df['classification_id'].unique()),
            height=15 * len(df['tax_id'].unique())
        ).transform_calculate(
            url=alt_kwargs.pop('url')
        ).mark_rect(
        ).encode(
            **alt_kwargs
        )

        if title:
            chart = chart.properties(title=title)

        if haxis:
            if return_chart:
                return top_label & chart
            else:
                (top_label & chart).display()
        else:
            if return_chart:
                return chart
            else:
                chart.interactive().display()
