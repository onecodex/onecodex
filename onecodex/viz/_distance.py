import numpy as np
import pandas as pd
import altair as alt
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from itertools import chain

from onecodex.exceptions import OneCodexException


class VizDistanceMixin():
    def plot_distance(self, rank='auto', metric='braycurtis',
                      label=None, title=None, xlabel=None, ylabel=None, tooltip=None):
        """Plot beta diversity distance matrix as a heatmap and dendrogram.

        metric -- one of: braycurtis, manhattan, jaccard, unifrac, unweighted_unifrac

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
            tooltip (list) -- display these metadata fields when points are hovered over
        """

        if rank is None:
            raise OneCodexException('Please specify a rank or \'auto\' to choose automatically')
        else:
            rank = self._get_auto_rank(rank)

        if len(self._results) < 2:
            raise OneCodexException('`plot_distance` requires 2 or more valid classification results.')

        # if taxonomy trees are inconsistent, unifrac will not work
        if metric in ('braycurtis', 'bray-curtis', 'bray curtis'):
            distances = self.beta_diversity(metric='braycurtis', rank=rank)
        elif metric in ('manhattan', 'cityblock'):
            distances = self.beta_diversity(metric='manhattan', rank=rank)
        elif metric == 'jaccard':
            distances = self.beta_diversity(metric='jaccard', rank=rank)
        elif metric in ('unifrac', 'weighted_unifrac'):
            distances = self.unifrac(weighted=True, rank=rank)
        elif metric == 'unweighted_unifrac':
            distances = self.unifrac(weighted=False, rank=rank)
        else:
            raise OneCodexException('Metric must be one of: braycurtis, manhattan, jaccard, '
                                    'weighted_unifrac, unweighted_unifrac')

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

        (dendro_chart | chart).display()
