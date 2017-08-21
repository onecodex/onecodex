import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from onecodex.exceptions import OneCodexException
from onecodex.helpers import normalize_classifications
from onecodex.distance import braycurtis, cityblock, jaccard, unifrac


def plot_distance(analyses, metric='braycurtis',
                  title=None, label=None, xlabel=None, ylabel=None,
                  field='readcount_w_children', rank='species'):
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

    normed_classifications, metadata = normalize_classifications(analyses, label=label)
    if len(normed_classifications) < 2:
        raise OneCodexException('`plot_distance` requires 2 or more valid classification results.')

    # there is no uniqueness constraint on metadata names
    # so plot by uuid, then replace the labels in the dataframe with their names
    uuids = {}
    sample_names = {}
    for idx, analysis in enumerate(normed_classifications):
        uuids[analysis.id] = analysis.id
        sample_names[analysis.id] = metadata.loc[idx, '_display_name']

    distances = f(normed_classifications, field=field, rank=rank)
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

    # Labels
    if xlabel is not None:
        plt.gca().set_xlabel(xlabel)
    if ylabel is not None:
        plt.gca().set_ylabel(ylabel)

    if title:
        g.fig.suptitle(title)
    plt.show()
