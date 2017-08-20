import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

from onecodex.helpers import collate_analysis_results, normalize_analyses


def plot_pca(analyses, title=None, threshold=None, field='readcount_w_children',
             rank=None, color_by=None, org_vectors=0, org_vectors_scale=None):
    # color_by: piece of metadata to color by
    # org_vectors: boolean; whether to plot the most highly contributing organisms
    # org_vectors_scale_factor: scale factor to modify the length of the organism vectors

    """Perform Principal Components Analysis to visualize the similarity of samples."""

    normed_analyses, metadata = normalize_analyses(analyses)
    df = collate_analysis_results(normed_analyses, field=field, rank=rank)

    # normalize the magnitude of the data
    df = (df.T / df.sum(axis=1)).T

    pca = PCA()
    pca_vals = pca.fit(df.values).transform(df.values)
    pca_vals = pd.DataFrame(pca_vals, index=df.index)
    pca_vals.rename(columns=lambda x: "PCA{}".format(x + 1), inplace=True)
    metadata.index = df.index
    plot_data = pd.concat([pca_vals, metadata.fillna('N/A')], axis=1)

    # Scatter plot of PCA
    sns.set(style="whitegrid")
    if color_by is None:
        g = sns.lmplot('PCA1', 'PCA2', data=plot_data, fit_reg=False)
    else:
        g = sns.lmplot('PCA1', 'PCA2', data=plot_data, fit_reg=False, hue=color_by)

    # Plot the organism eigenvectors that contribute the most
    if org_vectors > 0:
        # Plot the most highly contributing taxa
        magnitudes = np.sqrt(pca.components_[0] ** 2 + pca.components_[1] ** 2)
        magnitudes.sort()
        cutoff = magnitudes[-1 * org_vectors]
        # we can't use the "levels" method here b/c https://stackoverflow.com/questions/28772494
        tax_ids = [i[0] for i in df.columns.values]
        tax_id_map = {i[0]: i[1] for i in df.columns}
        if org_vectors_scale is None:
            org_vectors_scale = 0.8 * np.max(pca_vals.abs().values)
        for taxid, var1, var2 in \
                zip(tax_ids, pca.components_[0, :], pca.components_[1, :]):
            if np.sqrt(var1 ** 2 + var2 ** 2) >= cutoff:
                g.axes.flat[0].annotate("{} ({})".format(tax_id_map[taxid], taxid),
                                        xytext=(var1 * float(org_vectors_scale),
                                                var2 * float(org_vectors_scale)),
                                        xy=(0, 0), size=8,
                                        arrowprops={'facecolor': 'black',
                                                    'width': 1, 'headwidth': 5})
    plt.xlabel("PCA1 ({}%)".format(round(pca.explained_variance_ratio_[0] * 100, 2)))
    plt.ylabel("PCA2 ({}%)".format(round(pca.explained_variance_ratio_[1] * 100, 2)))
    if title:
        plt.title(title)
    plt.show()
