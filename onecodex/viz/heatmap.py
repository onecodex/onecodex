from onecodex.exceptions import OneCodexException
from onecodex.viz.helpers import collate_analysis_results, normalize_analyses


def plot_heatmap(analyses, title=None, top_n=20, threshold=None, field='readcount_w_children'):
    import matplotlib.pyplot as plt
    import seaborn as sns

    assert len(analyses) > 1

    if not (threshold or top_n):
        raise OneCodexException('Please set either "threshold" or "top_n"')

    normed_analyses, metadata = normalize_analyses(analyses)
    df = collate_analysis_results(normed_analyses, field=field)

    df.columns = ['{} ({})'.format(v[1], v[0]) for v in df.columns.values]
    df['display_name'] = [m.name if m.name is not None else m.sample.filename for m in metadata]
    df.set_index('display_name', inplace=True)

    if threshold and top_n:
        idx = df.sum(axis=0).sort_values(ascending=False).head(top_n).index
        df = df.loc[:, idx]  # can't use idx with multiple conditions
        g = sns.clustermap(df.loc[:, df.max() >= threshold])
    elif threshold:
        g = sns.clustermap(df.loc[:, df.max() >= threshold])
    elif top_n:
        idx = df.sum(axis=0).sort_values(ascending=False).head(top_n).index
        g = sns.clustermap(df.loc[:, idx])

    # Rotate the margin labels
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

    if title:
        g.fig.suptitle(title)

    plt.show()
