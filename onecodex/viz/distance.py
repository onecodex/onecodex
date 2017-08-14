import pandas as pd

from onecodex.viz.helpers import normalize_analyses, collate_analysis_results
from scipy.spatial.distance import braycurtis, cityblock
from onecodex.lib.diversity import unifrac, jaccard_dissimilarity


def plot_distance(analyses, title=None,
                  field='abundance', distance_metric='bray-curtis'):
    import matplotlib.pyplot as plt
    import seaborn as sns

    if distance_metric == 'bray-curtis':
        f = braycurtis
    elif distance_metric == 'manhattan':
        f = cityblock
    elif distance_metric == 'jaccard':
        f = jaccard_dissimilarity
    elif distance_metric == 'unifrac':
        f = unifrac
    else:
        print("Please specify 'distance_metric' as one "
              "of bray-curtis, manhattan, jaccard, or unifrac")
        return

    normed_analyses, metadata = normalize_analyses(analyses)
    df = collate_analysis_results(analyses, field=field)

    names = {}
    for analysis, idx in enumerate(normed_analyses):
        names[analysis.id] = metadata[idx]['name'] if metadata[idx]['name'] \
                                                   else metadata[idx]['filename']

    dists = {}
    for id1 in df.index.values:
        dists[names[id1]] = {}
        for id2 in df.index.values:
            dists[names[id1]][names[id2]] = f(df.loc[id1, :], df.loc[id2, :])
    dists = pd.DataFrame(dists)

    g = sns.clustermap(dists)
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    if title:
        g.fig.suptitle(title)
    plt.show()

    return dists
