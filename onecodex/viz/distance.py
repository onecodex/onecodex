import pandas as pd

from onecodex.exceptions import OneCodexException
from onecodex.viz.helpers import normalize_analyses
from onecodex.lib.distance_metrics import braycurtis, cityblock, unifrac, jaccard


def plot_distance(analyses, title=None, distance_metric='bray-curtis',
                  field='readcount_w_children', rank='species'):
    import matplotlib.pyplot as plt
    import seaborn as sns

    if distance_metric == 'bray-curtis':
        f = braycurtis
    elif distance_metric == 'manhattan':
        f = cityblock
    elif distance_metric == 'jaccard':
        f = jaccard
    elif distance_metric == 'unifrac':
        f = unifrac
    else:
        raise OneCodexException("Please specify 'distance_metric' as one "
                                "of bray-curtis, manhattan, jaccard, or unifrac")

    normed_analyses, metadata = normalize_analyses(analyses)

    names = {}
    for idx, analysis in enumerate(normed_analyses):
        # FIXME: if names are not unique, the dict will overwrite
        # names[analysis.id] = metadata[idx].name \
        #     if metadata[idx].name else analysis.sample.filename
        names[analysis.id] = analysis.id

    distances = f(analyses, field=field, rank=rank)
    ids = distances.ids
    distance_matrix = distances.data
    dists = {}
    for idx1, id1 in enumerate(ids):
        dists[names[id1]] = {}
        for idx2, id2 in enumerate(ids):
            dists[names[id1]][names[id2]] = distance_matrix[idx1][idx2]
    dists = pd.DataFrame(dists)

    g = sns.clustermap(dists)
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    if title:
        g.fig.suptitle(title)
    plt.show()

    return dists
