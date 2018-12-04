import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import altair as alt

from onecodex.exceptions import OneCodexException
from onecodex.helpers import collate_classification_results, normalize_classifications


def plot_pca(analyses,
             title=None, xlabel=None, ylabel=None, color=None, size=None, tooltip=None,
             field='readcount_w_children', rank=None, normalize=True):
    """Perform principal component analysis and plot first two axes.

    analyses (list) -- list of Samples, Classifications, or Analyses objects to be PCA'd

    field ('readcount_w_children' | 'readcount' | 'abundance')
        - 'readcount_w_children': total reads of this taxon and all its descendants
        - 'readcount': total reads of this taxon
        - 'abundance': genome size-normalized relative abundances, from shotgun sequencing
    rank (None | 'kingdom' | 'phylum' | 'class' | 'order' | 'family' | 'genus' | 'species')
        - None: include all ranks
        - 'kingdom' or others: restrict analysis to taxa at this rank
    normalize (bool): convert from read counts to relative abundances (each sample sums to 1.0)

    title (string) -- main title of the plot
    xlabel, ylabel (string) -- axes labels
    color (string) -- metadata field to color points by
    size (string) -- metadata field to size points by
        - For color and size, use 'taxid_N' where N is an arbitrary taxid to color or size
          points based on the abundance of that taxid
    tooltip (list) -- display these metadata fields when points are hovered over
        - Use 'taxid_N' where N is an arbitrary taxid to display its abundance in a tooltip
    """

    normed_classifications, metadata = normalize_classifications(analyses)
    df, tax_info = collate_classification_results(normed_classifications, field=field,
                                                  rank=rank, normalize=normalize)

    metadata.index = df.index
    metadata = metadata.fillna({field: 'N/A' for field in metadata.columns})

    if len(df) < 2:
        raise OneCodexException('`plot_pca` requires 2 or more valid classification results.')

    if tooltip:
        if not isinstance(tooltip, list):
            tooltip = [tooltip]
    else:
        tooltip = []

    # support point color/sizes and tooltips by taxid
    for param in [color, size] + tooltip:
        if param:
            if param.startswith('taxid_'):
                taxid = param[6:]

                if taxid not in df:
                    raise OneCodexException('Tax ID {} not found in analyses'.format(taxid))

                taxon_name = tax_info[taxid]['name']
                metadata[taxon_name] = df[taxid]

                if param == color:
                    color = taxon_name
                elif param == size:
                    size = taxon_name
                else:
                    tooltip[tooltip.index(param)] = taxon_name
            else:
                if param not in metadata and param != df.index.name:
                    raise OneCodexException('Column {} not found in metadata'.format(param))

    pca = PCA()
    pca_vals = pca.fit(df.values).transform(df.values)
    pca_vals = pd.DataFrame(pca_vals, index=df.index)
    pca_vals.rename(columns=lambda x: "PCA{}".format(x + 1), inplace=True)
    plot_data = pd.concat([pca_vals, metadata], axis=1).reset_index()

    # label the axes
    if xlabel is None:
        xlabel = 'PCA1 ({}%)'.format(round(pca.explained_variance_ratio_[0] * 100, 2))
    if ylabel is None:
        ylabel = 'PCA2 ({}%)'.format(round(pca.explained_variance_ratio_[1] * 100, 2))

    alt_kwargs = dict(
        x=alt.X('PCA1', axis=alt.Axis(title=xlabel)),
        y=alt.Y('PCA2', axis=alt.Axis(title=ylabel)),
        color=color,
        size=size,
        tooltip=tooltip,
        href='url:N',
        url='https://app.onecodex.com/classification/' + alt.datum.classification_id
    )

    # remove unused kwargs, altair doesn't like them
    for param in ('color', 'size', 'tooltip'):
        if not alt_kwargs[param]:
            alt_kwargs.pop(param)

    alt.renderers.enable('notebook')

    chart = alt.Chart(plot_data) \
               .transform_calculate(url=alt_kwargs.pop('url')) \
               .mark_circle() \
               .encode(**alt_kwargs)

    if title:
        chart = chart.properties(title=title)

    # Plot the organism eigenvectors that contribute the most
    if org_vectors > 0:
        plot_data = {
            'x': [],
            'y': [],
            'o': [],             # order these points should be connected in
            'Eigenvectors': []
        }

        # Plot the most highly contributing taxa
        magnitudes = np.sqrt(pca.components_[0] ** 2 + pca.components_[1] ** 2)
        magnitudes.sort()
        cutoff = magnitudes[-1 * org_vectors]

        if org_vectors_scale is None:
            org_vectors_scale = 0.8 * np.max(pca_vals.abs().values)

        for taxid, var1, var2 in zip(df.columns.values, pca.components_[0, :], pca.components_[1, :]):
            if np.sqrt(var1 ** 2 + var2 ** 2) >= cutoff:
                plot_data['x'].extend([0, var1 * float(org_vectors_scale)])
                plot_data['y'].extend([0, var2 * float(org_vectors_scale)])
                plot_data['o'].extend([0, 1])
                plot_data['Eigenvectors'].extend([tax_info[taxid]['name']] * 2)

                org_vectors -= 1

                if org_vectors == 0:
                    break

        plot_data = pd.DataFrame(plot_data)

        vector_chart = alt.Chart(plot_data) \
                          .mark_line(point=False) \
                          .encode(x=alt.X('x', axis=None),
                                  y=alt.Y('y', axis=None),
                                  order='o',
                                  color='Eigenvectors')

        (chart + vector_chart).interactive().display()
    else:
        chart.interactive().display()
