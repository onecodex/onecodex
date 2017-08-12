import matplotlib.pyplot as plt
import seaborn
import pandas
from sklearn.decomposition import PCA

from onecodex.viz.helpers import collate_analysis_results, normalize_analyses


def plot_pca(analyses, title=None, threshold=None, metric='abundance', color_by=None,
             org_vectors=False, org_vectors_scale_factor=8):
    # metric: 'abundance' or 'readcount'
    # color_by: piece of metadata to color by
    # org_vectors: boolean; whether to plot the most highly contributing organisms
    # org_vectors_scale_factor: scale factor to modify the length of the organism vectors

    """Perform Principal Components Analysis to visualize the similarity of samples."""

    normed_analyses, metadata = normalize_analyses(analyses)
    df = collate_analysis_results(normed_analyses, metric=metric)

    pca = PCA()
    pca_vals = pca.fit(df.values).transform(df.values)
    pca_vals = pandas.DataFrame(pca_vals, index=df.index)
    pca_vals.rename(columns=lambda x: "PCA{}".format(x + 1), inplace=True)

    pca = PCA()
    pca_vals = pca.fit(df.values).transform(df.values)
    pca_vals = pandas.DataFrame(pca_vals, index=df.index)
    pca_vals.rename(columns=lambda x: "PCA{}".format(x + 1), inplace=True)

    color_by_vals = []
    if color_by is not None:
        assert hasattr(metadata[0], color_by)
        # if the metadata does not have the color_by attribute, make it gray
        for md in metadata:
            if getattr(md, color_by) is None:
                color_by_vals.append('None')
            else:
                color_by_vals.append(getattr(md, color_by))

    pca_vals[color_by] = color_by_vals

    # Scatter plot of PCA
    if color_by is None:
        g = seaborn.lmplot('PCA1', 'PCA2', data=pca_vals, fit_reg=False)
    else:
        g = seaborn.lmplot('PCA1', 'PCA2', data=pca_vals, fit_reg=False, hue=color_by)

    # Plot the organism eigenvectors that contribute the most
    if org_vectors:
        # Plot the most highly contributing taxa
        n = 3  # Number of taxa to display
        cutoff = pandas.DataFrame(pca.components_[0:2, :]).abs().sum().sort_values().values[-1 * n]
        for taxid, var1, var2 in zip(df.columns, pca.components_[0, :], pca.components_[1, :]):
            if abs(var1) + abs(var2) >= cutoff:
                g.axes.flat[0].annotate("{} ({})".format(df['tax_name'], taxid),
                                        xytext=(var1/float(org_vectors_scale_factor),
                                                var2/float(org_vectors_scale_factor)),
                                        xy=(0, 0), size=8,
                                        arrowprops={'facecolor': 'black',
                                                    'width': 1, 'headwidth': 5})
    plt.xlabel("PCA1 ({}%)".format(round(pca.explained_variance_ratio_[0] * 100, 2)))
    plt.ylabel("PCA2 ({}%)".format(round(pca.explained_variance_ratio_[1] * 100, 2)))
    if title:
        plt.title(title)

    plt.show()
