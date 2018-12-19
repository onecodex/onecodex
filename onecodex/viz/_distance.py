import numpy as np
import pandas as pd
import altair as alt
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr
from skbio.stats import ordination
from sklearn import manifold
from sklearn.metrics.pairwise import euclidean_distances
from itertools import chain
import warnings

from onecodex.exceptions import OneCodexException
from onecodex.distance import DistanceMixin


class VizDistanceMixin(DistanceMixin):
    def _compute_distance(self, rank, metric):
        if rank is None:
            raise OneCodexException('Please specify a rank or \'auto\' to choose automatically')

        if len(self._results) < 2:
            raise OneCodexException('`plot_distance` requires 2 or more valid classification results.')

        # if taxonomy trees are inconsistent, unifrac will not work
        if callable(metric):
            distances = metric(self, rank=rank)
        elif metric in ('braycurtis', 'bray-curtis', 'bray curtis'):
            distances = self.beta_diversity(metric='braycurtis', rank=rank)
        elif metric in ('manhattan', 'cityblock'):
            distances = self.beta_diversity(metric='cityblock', rank=rank)
        elif metric == 'jaccard':
            distances = self.beta_diversity(metric='jaccard', rank=rank)
        elif metric in ('unifrac', 'weighted_unifrac'):
            distances = self.unifrac(weighted=True, rank=rank)
        elif metric == 'unweighted_unifrac':
            distances = self.unifrac(weighted=False, rank=rank)
        else:
            raise OneCodexException('Metric must be one of: braycurtis, manhattan, jaccard, '
                                    'weighted_unifrac, unweighted_unifrac')

        return distances

    def plot_distance(self, rank='auto', metric='braycurtis',
                      title=None, xlabel=None, ylabel=None, tooltip=None, return_chart=False):
        """Plot beta diversity distance matrix as a heatmap and dendrogram.

        Parameters
        ----------
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.
        metric : {'braycurtis', 'manhattan', 'jaccard', 'unifrac', 'unweighted_unifrac}, optional
            Function to use when calculating the distance between two samples.
        title : `string`, optional
            Text label at the top of the plot.
        xlabel : `string`, optional
            Text label along the horizontal axis.
        ylabel : `string`, optional
            Text label along the vertical axis.
        tooltip : `string` or `list`, optional
            A string or list containing strings representing metadata fields. When a point in the
            plot is hovered over, the value of the metadata associated with that sample will be
            displayed in a modal.

        Examples
        --------
        Plot the weighted UniFrac distance between all our samples, using counts at the genus level.

        >>> plot_distance(rank='genus', metric='unifrac')
        """
        distances = self._compute_distance(rank, metric)

        # this will be passed to the heatmap chart as a dataframe eventually
        plot_data = {
            '1) Label': [],
            '2) Label': [],
            'Distance': [],
            'classification_id': []
        }

        # here we figure out what to put in the tooltips and get the appropriate data
        if tooltip:
            if not isinstance(tooltip, list):
                tooltip = [tooltip]
        else:
            tooltip = []

        magic_metadata, magic_fields = self.magic_metadata_fetch(tooltip)
        formatted_fields = []

        for tip, magic_field in magic_fields.items():
            field_group = []

            for i in (1, 2):
                field = '{}) {}'.format(i, magic_field)
                plot_data[field] = []
                field_group.append(field)

            formatted_fields.append(field_group)

        # must convert to long format for heatmap plotting
        for idx1, id1 in enumerate(distances.ids):
            for idx2, id2 in enumerate(distances.ids):
                if idx1 == idx2:
                    plot_data['Distance'].append(np.nan)
                else:
                    plot_data['Distance'].append(distances.data[idx1][idx2])

                plot_data['1) Label'].append(self._metadata['_display_name'][id1])
                plot_data['2) Label'].append(self._metadata['_display_name'][id2])
                plot_data['classification_id'].append(id1)

                for field_group, magic_field in zip(formatted_fields, magic_fields.values()):
                    plot_data[field_group[0]].append(magic_metadata[magic_field][id1])
                    plot_data[field_group[1]].append(magic_metadata[magic_field][id2])

        plot_data = pd.DataFrame(data=plot_data)

        # turn the distances returned by skbio into a simple distance matrix
        dists = distances.to_data_frame()

        # here we use scipy to perform average-linkage clustering on the distance matrix
        clustering = hierarchy.linkage(squareform(dists), method='average')
        tree = hierarchy.dendrogram(clustering, no_plot=True)
        class_ids_in_order = [dists.index[int(x)] for x in tree['ivl']]
        names_in_order = self._metadata['_display_name'][class_ids_in_order].tolist()

        # it's important to tell altair to order the cells in the heatmap according to the clustering
        # obtained from scipy
        alt_kwargs = dict(
            x=alt.X('1) Label:N', axis=alt.Axis(title=xlabel), sort=names_in_order),
            y=alt.Y('2) Label:N', axis=alt.Axis(title=ylabel, orient='right'), sort=names_in_order),
            color='Distance:Q',
            tooltip=['1) Label', '2) Label', 'Distance:Q'] + list(chain.from_iterable(formatted_fields)),
            href='url:N',
            url='https://app.onecodex.com/classification/' + alt.datum.classification_id
        )

        chart = alt.Chart(plot_data,
                          width=15 * len(distances.ids),
                          height=15 * len(distances.ids)) \
                   .transform_calculate(url=alt_kwargs.pop('url')) \
                   .mark_rect() \
                   .encode(**alt_kwargs)

        if title:
            chart = chart.properties(title=title)

        # here we convert the dendrogram generated by scipy into a tree by plotting lines
        plot_data = {
            'x': [],
            'y': [],
            'o': [],  # order these points should be connected in
            'b': []   # one number per branch
        }

        for idx, (i, d) in enumerate(zip(tree['icoord'], tree['dcoord'])):
            plot_data['x'].extend(map(lambda x: -x, d))
            plot_data['y'].extend(map(lambda x: -x, i))
            plot_data['o'].extend([0, 1, 2, 3])
            plot_data['b'].extend([idx] * 4)

        plot_data = pd.DataFrame(plot_data)

        dendro_chart = alt.Chart(plot_data,
                                 width=100,
                                 height=15 * len(distances.ids)) \
                          .mark_line(point=False, opacity=0.5) \
                          .encode(x=alt.X('x', axis=None),
                                  y=alt.Y('y', axis=None),
                                  order='o',
                                  color=alt.Color('b:N',
                                                  scale=alt.Scale(domain=list(range(100)),
                                                                  range=['black'] * 100),
                                                  legend=None))

        if return_chart:
            return (dendro_chart | chart)
        else:
            (dendro_chart | chart).display()

    def plot_mds(self, rank='auto', metric='braycurtis', method='pcoa',
                 title=None, xlabel=None, ylabel=None, color=None, size=None, tooltip=None,
                 return_chart=False):
        """Plot beta diversity distance matrix using multidimensional scaling (MDS).

        Parameters
        ----------
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.
        metric : {'braycurtis', 'manhattan', 'jaccard', 'unifrac', 'unweighted_unifrac}, optional
            Function to use when calculating the distance between two samples.
        method : {'pcoa', 'smacof'}
            Algorithm to use for ordination. PCoA uses eigenvalue decomposition and is not well
            suited to non-euclidean distance functions. SMACOF is an iterative optimization strategy
            that can be used as an alternative.
        title : `string`, optional
            Text label at the top of the plot.
        xlabel : `string`, optional
            Text label along the horizontal axis.
        ylabel : `string`, optional
            Text label along the vertical axis.
        size : `string` or `tuple`, optional
            A string or a tuple containing strings representing metadata fields. The size of points
            in the resulting plot will change based on the metadata associated with each sample.
        color : `string` or `tuple`, optional
            A string or a tuple containing strings representing metadata fields. The color of points
            in the resulting plot will change based on the metadata associated with each sample.
        tooltip : `string` or `list`, optional
            A string or list containing strings representing metadata fields. When a point in the
            plot is hovered over, the value of the metadata associated with that sample will be
            displayed in a modal.

        Examples
        --------
        Scatter plot of weighted UniFrac distance between all our samples, using counts at the genus
        level.

        >>> plot_mds(rank='genus', metric='unifrac')

        Notes
        -----
        **For `smacof`**: The values reported on the axis labels are Pearson's correlations between
        the distances between points on each axis alone, and the corresponding distances in the
        distance matrix calculated using the user-specified metric. These values are related to the
        effectiveness of the MDS algorithm in placing points on the scatter plot in such a way that
        they truly represent the calculated distances. They do not reflect how well the distance
        metric captures similarities between the underlying data (in this case, an OTU table).
        """
        dists = self._compute_distance(rank, metric).to_data_frame()

        # here we figure out what to put in the tooltips and get the appropriate data
        if tooltip:
            if not isinstance(tooltip, list):
                tooltip = [tooltip]
        else:
            tooltip = []

        tooltip = list(set(['Label', color, size] + tooltip))

        magic_metadata, magic_fields = self.magic_metadata_fetch(tooltip)

        if method == 'smacof':
            # adapted from https://scikit-learn.org/stable/auto_examples/manifold/plot_mds.html
            x_field = 'MDS1'
            y_field = 'MDS2'

            seed = np.random.RandomState(seed=3)
            mds = manifold.MDS(
                max_iter=3000, eps=1e-12, random_state=seed, dissimilarity="precomputed", n_jobs=1
            )
            pos = mds.fit(dists).embedding_
            plot_data = pd.DataFrame(pos, columns=[x_field, y_field], index=dists.index)
            plot_data = plot_data.div(plot_data.abs().max(axis=0), axis=1)  # normalize to [0,1]

            # determine how much of the original distance is captured by each of the axes after MDS.
            # this implementation of MDS does not use eigen decomposition and so there's no simple
            # way of returning a 'percent of variance explained' value
            r_squared = []

            for axis in [0, 1]:
                mds_dist = pos.copy()
                mds_dist[::, axis] = 0
                mds_dist = squareform(euclidean_distances(mds_dist).round(6))
                r_squared.append(pearsonr(mds_dist, squareform(dists))[0])

            # label the axes
            x_extra_label = 'r² = %.02f' % (r_squared[0],)
            y_extra_label = 'r² = %.02f' % (r_squared[1],)
        elif method == 'pcoa':
            # suppress eigenvalue warning from skbio--not because it's an invalid warning, but
            # because lots of folks in the field run pcoa on these distances functions, even if
            # statistically inappropriate. perhaps this will change if we ever become more
            # opinionated about the analyses that we allow our users to do (roo)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ord_result = ordination.pcoa(dists.round(6))  # round to avoid float precision errors

            plot_data = ord_result.samples.iloc[:, [0, 1]]  # get first two components
            plot_data = plot_data.div(plot_data.abs().max(axis=0), axis=1)  # normalize to [0,1]
            plot_data.index = dists.index
            x_field, y_field = plot_data.columns.tolist()  # name of first two components

            x_extra_label = '%0.02f%%' % (ord_result.proportion_explained[0] * 100,)
            y_extra_label = '%0.02f%%' % (ord_result.proportion_explained[1] * 100,)
        else:
            raise OneCodexException('MDS method must be one of: smacof, pcoa')

        # label the axes
        if xlabel is None:
            xlabel = '{} ({})'.format(x_field, x_extra_label)
        if ylabel is None:
            ylabel = '{} ({})'.format(y_field, y_extra_label)

        plot_data = pd.concat([plot_data, magic_metadata], axis=1).reset_index()

        alt_kwargs = dict(
            x=alt.X(x_field, axis=alt.Axis(title=xlabel)),
            y=alt.Y(y_field, axis=alt.Axis(title=ylabel)),
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

        if return_chart:
            return chart
        else:
            chart.interactive().display()
