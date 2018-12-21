import altair as alt
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from sklearn.metrics.pairwise import euclidean_distances

from onecodex.exceptions import OneCodexException


class VizHeatmapMixin():
    def plot_heatmap(self, rank='auto', normalize='auto', top_n=20, threshold=None,
                     title=None, xlabel=None, ylabel=None, tooltip=None, return_chart=False,
                     linkage='average'):
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

        df = self.results(
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

        magic_metadata, magic_fields = self.magic_metadata_fetch(tooltip)

        # add columns for prettier display
        df['display_name'] = self._metadata['_display_name'][df['classification_id']].tolist()
        df['tax_name'] = ['{} ({})'.format(self._taxonomy['name'][t], t) for t in df['tax_id']]

        # and for metadata
        for f in tooltip:
            df[magic_fields[f]] = magic_metadata[magic_fields[f]][df['classification_id']].tolist()

        # use scipy to perform average-linkage clustering on euclidean distances (by taxa)
        df_for_clustering = self.results(
            rank=rank,
            normalize=normalize,
            top_n=top_n,
            threshold=threshold
        ).T
        taxa_dist = euclidean_distances(df_for_clustering).round(6)
        clustering = hierarchy.linkage(squareform(taxa_dist), method=linkage)
        tree = hierarchy.dendrogram(clustering, no_plot=True)
        tax_ids_in_order = [df_for_clustering.index[int(x)] for x in tree['ivl']]
        tax_names_in_order = ['{} ({})'.format(self._taxonomy['name'][t], t) for t in tax_ids_in_order]

        # use scipy to perform average-linkage clustering on euclidean distances (by sample)
        df_for_clustering = df_for_clustering.T
        sample_dist = euclidean_distances(df_for_clustering).round(6)
        clustering = hierarchy.linkage(squareform(sample_dist), method=linkage)
        tree = hierarchy.dendrogram(clustering, no_plot=True)
        c_ids_in_order = [df_for_clustering.index[int(x)] for x in tree['ivl']]
        labels_in_order = [self._metadata['_display_name'][t] for t in c_ids_in_order]

        alt_kwargs = dict(
            x=alt.X('display_name:N', axis=alt.Axis(title=xlabel), sort=labels_in_order),
            y=alt.Y('tax_name:N', axis=alt.Axis(title=ylabel), sort=tax_names_in_order),
            color='{}:Q'.format(self.field),
            tooltip=['{}:Q'.format(self.field)] + tooltip,
            href='url:N',
            url='https://app.onecodex.com/classification/' + alt.datum.classification_id
        )

        chart = alt.Chart(df,
                          width=15 * len(df['classification_id'].unique()),
                          height=15 * len(df['tax_id'].unique())) \
                   .transform_calculate(url=alt_kwargs.pop('url')) \
                   .mark_rect() \
                   .encode(**alt_kwargs)

        if title:
            chart = chart.properties(title=title)

        if return_chart:
            return chart
        else:
            chart.interactive().display()
