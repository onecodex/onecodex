import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import altair as alt

from onecodex.exceptions import OneCodexException


class VizPCAMixin():
    def plot_pca(self, rank='auto', normalize='auto', org_vectors=0, org_vectors_scale=None,
                 title=None, xlabel=None, ylabel=None, color=None, size=None, tooltip=None):
        """Perform principal component analysis and plot first two axes.

        Options for tabulation of classification results:
            rank ('auto' | kingdom' | 'phylum' | 'class' | 'order' | 'family' | 'genus' | 'species')
                - 'auto': choose automatically based on fields
                - 'kingdom' or others: restrict analysis to taxa at this rank
            normalize ('auto' | True | False):
                - 'auto': normalize data if readcount or readcount_w_children
                -  True: convert from read counts to relative abundances (each sample sums to 1.0)

        Options for plotting:
            title (string) -- main title of the plot
            xlabel, ylabel (string) -- axes labels
            size, color (string) -- metadata field to size or color points by
            tooltip (list) -- display these metadata fields when points are hovered over
                - For size, color, and tooltip: you can specify a tax_id or tax_name and its abundance
                  will be automatically inserted into the size/color/tooltip field
        """

        if rank is None:
            raise OneCodexException('Please specify a rank or \'auto\' to choose automatically')
        else:
            rank = self._get_auto_rank(rank)

        if len(self._results) < 2:
            raise OneCodexException('`plot_pca` requires 2 or more valid classification results.')

        if len(self._taxonomy) < 2:
            raise OneCodexException('Too few taxa in results. Need at least 2 for PCA.')

        df = self.results(rank=rank, normalize=normalize)

        if tooltip:
            if not isinstance(tooltip, list):
                tooltip = [tooltip]
        else:
            tooltip = []

        tooltip = list(set(['Label', color, size] + tooltip))

        magic_metadata, magic_fields = self.magic_metadata_fetch(tooltip)

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
        plot_data = pd.concat([pca_vals.ix[:, ('PCA1', 'PCA2')], magic_metadata], axis=1).reset_index()

        alt_kwargs = dict(
            x=alt.X('PCA1', axis=alt.Axis(title=xlabel)),
            y=alt.Y('PCA2', axis=alt.Axis(title=ylabel)),
            tooltip=[magic_fields[t] for t in tooltip if t],
            href='url:N',
            url='https://app.onecodex.com/classification/' + alt.datum.classification_id
        )

        # only add these parameters if they are in use
        if color:
            alt_kwargs['color'] = magic_fields[color]
        if size:
            alt_kwargs['size'] = magic_fields[size]

        chart = alt.Chart(plot_data) \
                   .transform_calculate(url=alt_kwargs.pop('url')) \
                   .mark_circle() \
                   .encode(**alt_kwargs)

        if title:
            chart = chart.properties(title=title)

        # plot the organism eigenvectors that contribute the most
        if org_vectors > 0:
            plot_data = {
                'x': [],
                'y': [],
                'o': [],             # order these points should be connected in
                'Eigenvectors': []
            }

            magnitudes = np.sqrt(pca.components_[0] ** 2 + pca.components_[1] ** 2)
            magnitudes.sort()
            cutoff = magnitudes[-1 * org_vectors]

            if org_vectors_scale is None:
                org_vectors_scale = 0.8 * np.max(pca_vals.abs().values)

            for tax_id, var1, var2 in zip(df.columns.values, pca.components_[0, :], pca.components_[1, :]):
                if np.sqrt(var1 ** 2 + var2 ** 2) >= cutoff:
                    plot_data['x'].extend([0, var1 * float(org_vectors_scale)])
                    plot_data['y'].extend([0, var2 * float(org_vectors_scale)])
                    plot_data['o'].extend([0, 1])
                    plot_data['Eigenvectors'].extend([self._taxonomy['name'][tax_id]] * 2)

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
