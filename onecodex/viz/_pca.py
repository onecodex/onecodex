import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import altair as alt
import warnings

from onecodex.exceptions import OneCodexException
from onecodex.helpers import collate_classification_results, normalize_classifications


def plot_pca(analyses, threshold=None,
             title=None, hue=None, xlabel=None, ylabel=None,
             org_vectors=0, org_vectors_scale=None,
             field='readcount_w_children', rank='genus', normalize=True,
             label=None, color=None, size=None, tooltip=None):
    """Perform principal component analysis and plot first two axes.

    analyses (list) -- list of Samples, Classifications, or Analyses objects to be PCA'd

    Options for tabulation of classification results:
        field ('readcount_w_children' | 'readcount' | 'abundance')
            - 'readcount_w_children': total reads of this taxon and all its descendants
            - 'readcount': total reads of this taxon
            - 'abundance': genome size-normalized relative abundances, from shotgun sequencing
        rank ('kingdom' | 'phylum' | 'class' | 'order' | 'family' | 'genus' | 'species')
            - 'kingdom' or others: restrict analysis to taxa at this rank
        normalize (bool): convert from read counts to relative abundances (each sample sums to 1.0)

    Options for plotting:
        label (string) -- metadata field to label samples with
        title (string) -- main title of the plot
        xlabel, ylabel (string) -- axes labels
        size, color (string) -- metadata field to size or color points by
        tooltip (list) -- display these metadata fields when points are hovered over
            - For size, color, and tooltip: you can specify a tax_id or tax_name and its abundance
              will be automatically inserted into the size/color/tooltip field
    """

    # TODO: remove this code to rename hue to color
    if hue:
        color = hue

    # TODO: remove this code when org_vectors kwargs are removed
    if org_vectors or org_vectors_scale:
        warnings.warn('`plot_pca` no longer supports organism vectors')

    if rank is None:
        raise OneCodexException('Please specify a taxonomic rank')

    if not isinstance(analyses, list) or len(analyses) < 2:
        raise OneCodexException('`plot_pca` requires 2 or more valid classification results.')

    normed_classifications, metadata = normalize_classifications(analyses, label=label)
    df, tax_info = collate_classification_results(normed_classifications, field=field,
                                                  rank=rank, normalize=normalize)

    if len(tax_info) < 2:
        raise OneCodexException(
            'Two few cols (taxa) in classification results provided. Need at least 2 for PCA.'
        )

    metadata.index = df.index
    metadata = metadata.fillna({field: 'N/A' for field in metadata.columns})

    if tooltip:
        if not isinstance(tooltip, list):
            tooltip = [tooltip]
    else:
        tooltip = []

    # support point color/sizes and tooltips by taxid
    renamed_fields = {}

    for param in [color, size] + tooltip:
        if param:
            if param in metadata or param == df.index.name:
                # field is simply in the metadata dataframe
                renamed_fields[param] = param
                continue
            elif str(param) in df.keys():
                # it's a tax_id
                new_field_name = '{} ({})'.format(tax_info[str(param)]['name'], param)
                metadata[new_field_name] = df[str(param)]
                renamed_fields[param] = new_field_name
            else:
                # it might be a tax_name
                hits = []

                for t in tax_info:
                    if param.lower() in tax_info[t]['name'].lower():
                        hits.append(t)

                # take the lowest tax_id that maches the search query
                hits = sorted(hits, key=int)

                if hits:
                    new_field_name = '{} ({})'.format(tax_info[hits[0]]['name'], hits[0])
                    metadata[new_field_name] = df[hits[0]]
                    renamed_fields[param] = new_field_name
                else:
                    raise OneCodexException('Column {} not found in metadata'.format(param))

    pca = PCA()
    pca_vals = pca.fit(df.values).transform(df.values)
    pca_vals = pd.DataFrame(pca_vals, index=df.index)
    pca_vals.rename(columns=lambda x: "PCA{}".format(x + 1), inplace=True)

    # label the axes
    if xlabel is None:
        xlabel = 'PCA1 ({}%)'.format(round(pca.explained_variance_ratio_[0] * 100, 2))
    if ylabel is None:
        ylabel = 'PCA2 ({}%)'.format(round(pca.explained_variance_ratio_[1] * 100, 2))

    # don't send all the data to vega, just what we're plotting
    plot_data = pca_vals.ix[:, ('PCA1', 'PCA2')]

    for param in [color, size] + tooltip:
        if param and param != df.index.name:
            plot_data[renamed_fields[param]] = metadata[renamed_fields[param]]

    plot_data['Label'] = metadata['_display_name']
    plot_data = plot_data.reset_index()

    alt_kwargs = dict(
        x=alt.X('PCA1', axis=alt.Axis(title=xlabel)),
        y=alt.Y('PCA2', axis=alt.Axis(title=ylabel)),
        color=renamed_fields[color] if color else None,
        size=renamed_fields[size] if size else None,
        # put the sample label in tooltip by default
        tooltip=['Label'] + [renamed_fields[t] for t in tooltip if t],
        href='url:N',
        url='https://app.onecodex.com/classification/' + alt.datum.classification_id
    )

    # remove unused kwargs, altair doesn't like them
    for param in ['color', 'size', 'tooltip']:
        if not alt_kwargs[param]:
            alt_kwargs.pop(param)

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
