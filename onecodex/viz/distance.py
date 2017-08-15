import pandas as pd

from onecodex.exceptions import OneCodexException
from onecodex.viz.helpers import normalize_analyses
from onecodex.lib.distance_metrics import braycurtis, cityblock, jaccard


def plot_distance(analyses, title=None, distance_metric='bray-curtis',
                  field='readcount_w_children', rank='species'):
    import matplotlib.pyplot as plt
    import seaborn as sns

    # if taxonomy trees are inconsistent, unifrac will not work
    if distance_metric == 'bray-curtis':
        f = braycurtis
    elif distance_metric in ['manhattan', 'cityblock']:
        f = cityblock
    elif distance_metric == 'jaccard':
        f = jaccard
    else:
        raise OneCodexException("Please specify 'distance_metric' as one "
                                "of bray-curtis, manhattan, or jaccard")

    normed_analyses, metadata = normalize_analyses(analyses)

    # there is no uniqueness constraint on metadata names
    # so plot by uuid, then replace the labels in the dataframe with their names
    uuids = {}
    sample_names = {}
    for idx, analysis in enumerate(normed_analyses):
        uuids[analysis.id] = analysis.id
        sample_names[analysis.id] = metadata[idx].name \
            if metadata[idx].name else analysis.sample.filename

    distances = f(analyses, field=field, rank=rank)
    ids = distances.ids
    distance_matrix = distances.data
    dists = {}
    for idx1, id1 in enumerate(ids):
        dists[uuids[id1]] = {}
        for idx2, id2 in enumerate(ids):
            dists[uuids[id1]][uuids[id2]] = distance_matrix[idx1][idx2]
    dists = pd.DataFrame(dists).rename(index=sample_names, columns=sample_names)
    g = sns.clustermap(dists)
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    if title:
        g.fig.suptitle(title)
    plt.show()

    return dists
