import altair as alt

from onecodex.exceptions import OneCodexException


class VizHeatmapMixin():
    def plot_heatmap(self, rank='auto', normalize='auto',
                     top_n=20, threshold=None,
                     title=None, xlabel=None, ylabel=None, tooltip=None):
        """Plot heatmap of taxa abundance/count data for several samples.

        Specify at least one of the following options:
            top_n (int) -- display the top N most abundant taxa
            threshold (float) -- display taxa more abundant than this threshold in 1 or more sample

        Options for tabulation of classification results:
            rank ('auto' | kingdom' | 'phylum' | 'class' | 'order' | 'family' | 'genus' | 'species')
                - 'auto': choose automatically based on fields
                - 'kingdom' or others: restrict analysis to taxa at this rank
            normalize ('auto' | True | False):
                - 'auto': normalize data if readcount or readcount_w_children
                -  True: convert from read counts to relative abundances (each sample sums to 1.0)

        Options for plotting:
            title (string) -- main title of the plot
            xlabel, ylabel (string) -- axes labels
            tooltip (list) -- display these metadata fields when points are hovered over
        """

        if rank is None:
            raise OneCodexException('Please specify a rank or \'auto\' to choose automatically')
        else:
            rank = self._get_auto_rank(rank)

        if not (threshold or top_n):
            raise OneCodexException('Please specify at least one of: threshold, top_n')

        if len(self._results) < 2:
            raise OneCodexException('`plot_heatmap` requires 2 or more valid classification results.')

        if rank == 'auto':
            if self.field == 'abundance':
                rank = 'species'
            else:
                rank = 'genus'

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

        alt_kwargs = dict(
            x=alt.X('display_name:N', axis=alt.Axis(title=xlabel)),
            y=alt.Y('tax_name:N', axis=alt.Axis(title=ylabel)),
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

        chart.interactive().display()
