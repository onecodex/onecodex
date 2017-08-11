import matplotlib as plt
import seaborn

from onecodex.viz.helpers import collate_analysis_results, normalize_analyses


def plot_heatmap(analyses, title=None, min_abund=0.01, metric='abundance'):
    # min_abund: Only plot taxa that reach this abundance in at least one sample
    # metric: 'abundance' or 'readcount'

    a, metadata = normalize_analyses(analyses)
    df = collate_analysis_results(a, metric=metric)

    # Only keep the metadata for samples that have microbial abundance data
    # metadata = metadata.loc[[ix for ix in df.index.values if ix in metadata.index.values], :]
    # metadata = [ix for ix in df.index.values if ix in metadata.index.values]

    # Order the abundance data the same as the metadata
    # df = df.loc[metadata.index.values, :]

    # Name the taxa
    df.rename(columns=lambda t: "{} ({})".format(df['tax_name'], t), inplace=True)

    # Name the samples (show the name if specified, otherwise fall back to the filename)
    df['display_name'] = [r['name'] if r['name'] is not None else r['filename']
                          for ix, r in metadata.iterrows()]
    df.set_index('display_name', inplace=True)

    # Plot a heatmap of the raw abundances
    g = seaborn.clustermap(df.loc[:, df.max() >= min_abund])

    # Rotate the margin labels
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

    if title:
        g.fig.suptitle(title)

    plt.show()
