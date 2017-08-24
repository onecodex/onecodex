import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy
import warnings

from onecodex.exceptions import OneCodexException
from onecodex.helpers import normalize_classifications
from onecodex.distance import braycurtis, cityblock, jaccard, unifrac


def plot_distance(analyses, metric='braycurtis',
                  title=None, label=None, xlabel=None, ylabel=None,
                  field='readcount_w_children', rank='species', **kwargs):
    """Plot beta diversity distance matrix.

    Additional **kwargs are passed to Seaborn's `sns.clustermap`.
    """
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

    sns.set(style=kwargs.pop('style', 'darkgrid'))

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

    # Plot cluster map; ignore new SciPy cluster warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', scipy.cluster.hierarchy.ClusterWarning)
        g = sns.clustermap(dists, **kwargs)

    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

    # Labels
    if xlabel is not None:
        plt.gca().set_xlabel(xlabel)
    if ylabel is not None:
        plt.gca().set_ylabel(ylabel)

    if title:
        g.fig.suptitle(title)
    plt.show()
