import pandas as pd
import altair as alt

from onecodex.exceptions import OneCodexException
from onecodex.helpers import collate_classification_results, normalize_classifications


def plot_heatmap(analyses, top_n=20, threshold=None,
                 title=None, label=None, xlabel=None, ylabel=None, tooltip=None,
                 field='readcount_w_children', rank='genus', normalize=True):
    """Plot heatmap of taxa abundance/count data for several samples.

    analyses (list) -- list of Samples, Classifications, or Analyses objects to be plotted

    Specify only one of the following options:
        top_n (int) -- display the top N most abundant taxa
        threshold (float) -- display only taxa more abundant than this threshold

    Options for tabulation of classification results:
        field ('readcount_w_children' | 'readcount' | 'abundance')
            - 'readcount_w_children': total reads of this taxon and all its descendants
            - 'readcount': total reads of this taxon
            - 'abundance': genome size-normalized relative abundances, from shotgun sequencing
        rank ('kingdom' | 'phylum' | 'class' | 'order' | 'family' | 'genus' | 'species')
            - 'kingdom' or others: restrict analysis to taxa at this rank
        normalize (bool): convert from read counts to relative abundances (each sample sums to 1.0)

    Options for plotting:
        label (string) -- metadata field to label samples with
        title (string) -- main title of the plot
        xlabel, ylabel (string) -- axes labels
        tooltip (list) -- display these metadata fields when points are hovered over
    """

    if rank is None:
        raise OneCodexException('Please specify a taxonomic rank')

    if not (threshold or top_n):
        raise OneCodexException('Please specify one of: threshold, top_n')

    if not isinstance(analyses, list) or len(analyses) < 2:
        raise OneCodexException('`plot_heatmap` requires 2 or more valid classification results.')

    normed_classifications, metadata = normalize_classifications(analyses, label=label)
    df, tax_info = collate_classification_results(normed_classifications, field=field,
                                                  rank=rank, normalize=normalize)

    metadata.index = df.index

    df.columns = ['{} ({})'.format(tax_info[tax_id]['name'], tax_id) for tax_id in df.columns.values]

    # filter out taxa to plot
    if top_n:
        idx = df.sum(axis=0).sort_values(ascending=False).head(top_n).index
        df = df.loc[:, idx]

    if threshold:
        df = df.loc[:, df.max() >= threshold]

    # prepare tooltips
    if tooltip:
        if not isinstance(tooltip, list):
            tooltip = [tooltip]
    else:
        tooltip = []

    for param in tooltip:
        if param not in metadata and param != df.index.name:
            raise OneCodexException('Column {} not found in metadata'.format(param))

    # transfer data into something altair can handle
    plot_data = {
        'classification_id': [],
        'display_name': [],
        'taxon': [],
        'value': []
    }

    for taxon in df:
        for sample, abundance in df[taxon].iteritems():
            plot_data['classification_id'].append(sample)
            plot_data['display_name'].append(metadata['_display_name'][sample])
            plot_data['taxon'].append(taxon)
            plot_data['value'].append(abundance)

            for tip in tooltip:
                try:
                    plot_data[tip].append(metadata[tip][sample])
                except KeyError:
                    plot_data[tip] = [metadata[tip][sample]]

    plot_data = pd.DataFrame(data=plot_data)

    alt_kwargs = dict(
        x=alt.X('display_name:N', axis=alt.Axis(title=xlabel)),
        y=alt.Y('taxon:N', axis=alt.Axis(title=ylabel)),
        color='value:Q',
        tooltip=['value:Q'] + tooltip,
        href='url:N',
        url='https://app.onecodex.com/classification/' + alt.datum.classification_id
    )

    chart = alt.Chart(plot_data,
                      width=15 * len(df.index),
                      height=15 * len(df.keys())) \
               .transform_calculate(url=alt_kwargs.pop('url')) \
               .mark_rect() \
               .encode(**alt_kwargs)

    if title:
        chart = chart.properties(title=title)

    chart.interactive().display()
