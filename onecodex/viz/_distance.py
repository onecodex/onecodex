import numpy as np
import pandas as pd
import altair as alt
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

from onecodex.exceptions import OneCodexException
from onecodex.helpers import normalize_classifications, collate_classification_results
from onecodex.distance import braycurtis, cityblock, jaccard, unifrac


def plot_distance(analyses, metric='braycurtis',
                  title=None, label=None, xlabel=None, ylabel=None,
                  field='readcount_w_children', rank='species', normalize=True):
    """Plot beta diversity distance matrix."""

    # if taxonomy trees are inconsistent, unifrac will not work
    if metric in ['braycurtis', 'bray-curtis', 'bray curtis']:
        f = braycurtis
    elif metric in ['manhattan', 'cityblock']:
        f = cityblock
    elif metric == 'jaccard':
        f = jaccard
    elif metric == 'unifrac':
        f = unifrac
    else:
        raise OneCodexException("'metric' must be one of "
                                "braycurtis, manhattan, jaccard, or unifrac")

    if rank is None:
        raise OneCodexException('Please specify a taxonomic rank')

    if not isinstance(analyses, list) or len(analyses) < 2:
        raise OneCodexException('`plot_distance` requires 2 or more valid classification results.')

    normed_classifications, metadata = normalize_classifications(analyses, label=label)
    df, tax_info = collate_classification_results(normed_classifications, field=field,
                                                  rank=rank, normalize=normalize)

    metadata.index = df.index
    distances = f(normed_classifications, field=field, rank=rank)

    plot_data = {
        'label1': [],
        'label2': [],
        'distance': []
    }

    dists = {}

    for idx1, id1 in enumerate(distances.ids):
        dists[id1] = {}

        for idx2, id2 in enumerate(distances.ids):
            if idx1 == idx2:
                plot_data['distance'].append(np.nan)
            else:
                plot_data['distance'].append(distances.data[idx1][idx2])

            plot_data['label1'].append(metadata['_display_name'][id1])
            plot_data['label2'].append(metadata['_display_name'][id2])

            dists[id1][id2] = distances.data[idx1][idx2]

    plot_data = pd.DataFrame(data=plot_data)

    dists = pd.DataFrame(dists)
    dists.index.name = 'classification_id'

    clustering = hierarchy.linkage(squareform(dists), method='average')
    tree = hierarchy.dendrogram(clustering, no_plot=True)
    class_ids_in_order = [dists.index[int(x)] for x in tree['ivl']]
    names_in_order = metadata['_display_name'][class_ids_in_order].tolist()

    alt_kwargs = dict(
        x=alt.X('label1:N', axis=alt.Axis(title=xlabel), sort=names_in_order),
        y=alt.Y('label2:N', axis=alt.Axis(title=ylabel, orient='right'), sort=names_in_order),
        color='distance:Q',
        tooltip=['label1', 'label2', 'distance:Q'],
    )

    chart = alt.Chart(plot_data,
                      width=15 * len(distances.ids),
                      height=15 * len(distances.ids)) \
               .mark_rect() \
               .encode(**alt_kwargs)

    if title:
        chart = chart.properties(title=title)

    plot_data = {
        'x': [],
        'y': [],
        'o': [],  # order these points should be connected in
        'b': []   # one number per branch
    }

    for idx, (i, d) in enumerate(zip(tree['icoord'], tree['dcoord'])):
        plot_data['x'].extend(map(lambda x: -x, d))
        plot_data['y'].extend(map(lambda x: -x, i))
        plot_data['o'].extend([0, 1, 2, 3])
        plot_data['b'].extend([idx] * 4)

    plot_data = pd.DataFrame(plot_data)

    dendro_chart = alt.Chart(plot_data,
                             width=100,
                             height=15 * len(distances.ids)) \
                      .mark_line(point=False, opacity=0.5) \
                      .encode(x=alt.X('x', axis=None),
                              y=alt.Y('y', axis=None),
                              order='o',
                              color=alt.Color('b:N',
                                              scale=alt.Scale(domain=list(range(100)),
                                                              range=['black'] * 100),
                                              legend=None))

    (dendro_chart | chart).display()
