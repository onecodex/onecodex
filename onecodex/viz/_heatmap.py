import matplotlib.pyplot as plt
import seaborn as sns

from onecodex.exceptions import OneCodexException
from onecodex.helpers import collate_classification_results, normalize_classifications


def plot_heatmap(analyses, top_n=20, threshold=None,
                 title=None, label=None, xlabel=None, ylabel=None,
                 field='readcount_w_children', rank=None, normalize=False, **kwargs):
    """Plot a heatmap of classification results by sample.

    Additional **kwargs are passed to Seaborn's `sns.clustermap`.
    """
    if not (threshold or top_n):
        raise OneCodexException('Please set either "threshold" or "top_n"')

    normed_classifications, metadata = normalize_classifications(analyses, label=label)
    df, tax_info = collate_classification_results(normed_classifications, field=field, rank=rank,
                                                  normalize=normalize)

    if len(df) < 2:
        raise OneCodexException('`plot_heatmap` requires 2 or more valid classification results.')

    df.columns = ['{} ({})'.format(tax_info[tax_id]['name'], tax_id) for tax_id in df.columns.values]
    df.index = metadata.loc[:, '_display_name']
    df.index.name = ''

    sns.set(style=kwargs.pop('style', 'darkgrid'))

    if threshold and top_n:
        idx = df.sum(axis=0).sort_values(ascending=False).head(top_n).index
        df = df.loc[:, idx]  # can't use idx with multiple conditions
        g = sns.clustermap(df.loc[:, df.max() >= threshold], **kwargs)
    elif threshold:
        g = sns.clustermap(df.loc[:, df.max() >= threshold], **kwargs)
    elif top_n:
        idx = df.sum(axis=0).sort_values(ascending=False).head(top_n).index
        g = sns.clustermap(df.loc[:, idx], **kwargs)

    # Rotate the margin labels
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

    # Labels
    if xlabel is not None:
        plt.gca().set_xlabel(xlabel)
    if ylabel is not None:
        plt.gca().set_ylabel(ylabel)

    if title:
        g.fig.suptitle(title)

    plt.show()
