import pandas as pd

from onecodex.exceptions import OneCodexException
from onecodex.viz.helpers import normalize_analyses, collate_analysis_results
from onecodex.lib.diversity import bray_curtis, cityblock, unifrac, jaccard_dissimilarity


def plot_distance(analyses, title=None, distance_metric='bray-curtis',
                  field='readcount_w_children', rank='species'):
    import matplotlib.pyplot as plt
    import seaborn as sns

    if distance_metric == 'bray-curtis':
        f = bray_curtis
    elif distance_metric == 'manhattan':
        f = cityblock
    elif distance_metric == 'jaccard':
        f = jaccard_dissimilarity
    elif distance_metric == 'unifrac':
        f = unifrac
    else:
        raise OneCodexException("Please specify 'distance_metric' as one "
                                "of bray-curtis, manhattan, jaccard, or unifrac")

    normed_analyses, metadata = normalize_analyses(analyses)
    df = collate_analysis_results(analyses, field=field)

    names = {}
    for analysis, idx in enumerate(normed_analyses):
        names[analysis.id] = metadata[idx]['name'] \
            if metadata[idx]['name'] else metadata[idx]['filename']

    dists = {}
    for id1 in df.index.values:
        dists[names[id1]] = {}
        for id2 in df.index.values:
            dists[names[id1]][names[id2]] = f(df.loc[id1, :], df.loc[id2, :],
                                              field=field, rank=rank)
    dists = pd.DataFrame(dists)

    g = sns.clustermap(dists)
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    if title:
        g.fig.suptitle(title)
    plt.show()

    return dists
