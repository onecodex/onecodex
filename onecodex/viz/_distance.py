# -*- coding: utf-8 -*-
from itertools import chain
import warnings

from onecodex.lib.enums import BetaDiversityMetric, Rank, Linkage, OrdinationMethod
from onecodex.exceptions import OneCodexException, PlottingException
from onecodex.distance import DistanceMixin
from onecodex.viz._primitives import (
    interleave_palette,
    prepare_props,
    get_classification_url,
    open_links_in_new_tab,
)
from onecodex.utils import is_continuous, has_missing_values


class VizDistanceMixin(DistanceMixin):
    def _compute_distance(self, rank, metric):
        if rank is None:
            raise OneCodexException("Please specify a rank or 'auto' to choose automatically")

        # if taxonomy trees are inconsistent, unifrac will not work
        if callable(metric):
            distances = metric(self, rank=rank)
        elif metric in (BetaDiversityMetric.BrayCurtis, "bray-curtis", "bray curtis"):
            distances = self.beta_diversity(metric=BetaDiversityMetric.BrayCurtis, rank=rank)
        elif metric in ("manhattan", BetaDiversityMetric.CityBlock):
            distances = self.beta_diversity(metric=BetaDiversityMetric.CityBlock, rank=rank)
        elif metric == BetaDiversityMetric.Jaccard:
            distances = self.beta_diversity(metric=BetaDiversityMetric.Jaccard, rank=rank)
        elif metric == BetaDiversityMetric.WeightedUnifrac:
            distances = self.unifrac(weighted=True, rank=rank)
        elif metric == BetaDiversityMetric.UnweightedUnifrac:
            distances = self.unifrac(weighted=False, rank=rank)
        elif metric == BetaDiversityMetric.Aitchison:
            distances = self.beta_diversity(metric=BetaDiversityMetric.Aitchison, rank=rank)
        else:
            raise OneCodexException(
                "Metric must be one of: {}".format(", ".join(BetaDiversityMetric.values()))
            )

        return distances

    def _cluster_by_sample(
        self, rank=Rank.Auto, metric=BetaDiversityMetric.BrayCurtis, linkage=Linkage.Average
    ):
        from scipy.cluster import hierarchy
        from scipy.spatial.distance import squareform
        from sklearn.metrics.pairwise import euclidean_distances

        if metric == "euclidean":
            dist_matrix = euclidean_distances(self._results).round(6)
        else:
            dist_matrix = self._compute_distance(rank=rank, metric=metric).to_data_frame().round(6)
        clustering = hierarchy.linkage(squareform(dist_matrix), method=linkage)
        scipy_tree = hierarchy.dendrogram(clustering, no_plot=True)
        ids_in_order = [self._results.index[int(x)] for x in scipy_tree["ivl"]]

        return {
            "dist_matrix": dist_matrix,
            "clustering": clustering,
            "scipy_tree": scipy_tree,
            "ids_in_order": ids_in_order,
        }

    def _cluster_by_taxa(self, linkage=Linkage.Average):
        from scipy.cluster import hierarchy
        from scipy.spatial.distance import squareform
        from sklearn.metrics.pairwise import euclidean_distances

        dist_matrix = euclidean_distances(self._results.T).round(6)
        clustering = hierarchy.linkage(squareform(dist_matrix), method=linkage)
        scipy_tree = hierarchy.dendrogram(clustering, no_plot=True)
        ids_in_order = [self._results.T.index[int(x)] for x in scipy_tree["ivl"]]
        labels_in_order = ["{} ({})".format(self.taxonomy["name"][t], t) for t in ids_in_order]

        return {
            "dist_matrix": dist_matrix,
            "clustering": clustering,
            "scipy_tree": scipy_tree,
            "ids_in_order": ids_in_order,
            "labels_in_order": labels_in_order,
        }

    def plot_distance(
        self,
        rank=Rank.Auto,
        metric=BetaDiversityMetric.BrayCurtis,
        title=None,
        xlabel=None,
        ylabel=None,
        tooltip=None,
        return_chart=False,
        linkage=Linkage.Average,
        label=None,
        width=None,
        height=None,
    ):
        """Plot beta diversity distance matrix as a heatmap and dendrogram.

        Parameters
        ----------
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.
        metric : {'braycurtis', 'cityblock', 'manhattan', 'jaccard', 'unifrac', 'unweighted_unifrac', 'aitchison'}, optional
            Function to use when calculating the distance between two samples.
            Note that 'cityblock' and 'manhattan' are equivalent metrics.
        linkage : {'average', 'single', 'complete', 'weighted', 'centroid', 'median'}
            The type of linkage to use when clustering axes.
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
        label : `string` or `callable`, optional
            A metadata field (or function) used to label each analysis. If passing a function, a
            dict containing the metadata for each analysis is passed as the first and only
            positional argument. The callable function must return a string.

        Examples
        --------
        Plot the weighted UniFrac distance between all our samples, using counts at the genus level.

        >>> plot_distance(rank='genus', metric='unifrac')
        """
        import altair as alt
        import numpy as np
        import pandas as pd
        from onecodex.viz import dendrogram

        if len(self._results) < 2:
            raise PlottingException(
                "There are too few samples for distance matrix plots after filtering. Please "
                "select 2 or more samples to plot."
            )

        # this will be passed to the heatmap chart as a dataframe eventually
        plot_data = {"1) Label": [], "2) Label": [], "Distance": [], "classification_id": []}

        # here we figure out what to put in the tooltips and get the appropriate data
        if tooltip:
            if not isinstance(tooltip, list):
                tooltip = [tooltip]
        else:
            tooltip = []

        tooltip.insert(0, "Label")

        magic_metadata, magic_fields = self._metadata_fetch(tooltip, label=label)
        formatted_fields = []

        for _, magic_field in magic_fields.items():
            field_group = []

            for i in (1, 2):
                field = "{}) {}".format(i, magic_field)
                plot_data[field] = []
                field_group.append(field)

            formatted_fields.append(field_group)

        clust = self._cluster_by_sample(rank=rank, metric=metric, linkage=linkage)

        # must convert to long format for heatmap plotting
        for idx1, id1 in enumerate(clust["dist_matrix"].index):
            for idx2, id2 in enumerate(clust["dist_matrix"].index):
                if idx1 == idx2:
                    plot_data["Distance"].append(np.nan)
                else:
                    plot_data["Distance"].append(clust["dist_matrix"].iloc[idx1, idx2])

                plot_data["classification_id"].append(id1)

                for field_group, magic_field in zip(formatted_fields, magic_fields.values()):
                    plot_data[field_group[0]].append(magic_metadata[magic_field][id1])
                    plot_data[field_group[1]].append(magic_metadata[magic_field][id2])

        plot_data = pd.DataFrame(data=plot_data)
        labels_in_order = magic_metadata["Label"][clust["ids_in_order"]].tolist()
        plot_data["url"] = plot_data["classification_id"].apply(get_classification_url)

        # it's important to tell altair to order the cells in the heatmap according to the clustering
        # obtained from scipy
        alt_kwargs = dict(
            x=alt.X("1) Label:N", axis=alt.Axis(title=xlabel), sort=labels_in_order),
            y=alt.Y(
                "2) Label:N", axis=alt.Axis(title=ylabel, orient="right"), sort=labels_in_order
            ),
            color=alt.Color("Distance:Q", legend=alt.Legend(title="Distance")),
            tooltip=list(chain.from_iterable(formatted_fields)) + ["Distance:Q"],
            href="url:N",
        )

        chart = (
            alt.Chart(
                plot_data,
                width=15 * len(clust["dist_matrix"].index),
                height=15 * len(clust["dist_matrix"].index),
            )
            .mark_rect()
            .encode(**alt_kwargs)
        )

        chart = chart.properties(**prepare_props(height=height, width=width))

        dendro_chart = dendrogram(clust["scipy_tree"])

        if height:
            cell_height = height / len(clust["dist_matrix"].index)
            dendro_chart = dendro_chart.properties(height=height - cell_height / 2)

        title_kwargs = prepare_props(title=title)
        concat_chart = alt.hconcat(dendro_chart, chart, spacing=0, **title_kwargs).configure_view(
            strokeWidth=0
        )
        open_links_in_new_tab(concat_chart)

        if return_chart:
            return concat_chart
        else:
            concat_chart.display()

    def plot_pcoa(self, *args, **kwargs):
        return self.plot_mds(*args, method=OrdinationMethod.Pcoa, **kwargs)

    def plot_mds(
        self,
        rank=Rank.Auto,
        metric=BetaDiversityMetric.BrayCurtis,
        method=OrdinationMethod.Pcoa,
        title=None,
        xlabel=None,
        ylabel=None,
        color=None,
        size=None,
        tooltip=None,
        return_chart=False,
        label=None,
        mark_size=100,
        width=None,
        height=None,
    ):
        """Plot beta diversity distance matrix using multidimensional scaling (MDS).

        Parameters
        ----------
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.
        metric : {'braycurtis', 'cityblock', 'manhattan', 'jaccard', 'unifrac', 'unweighted_unifrac', 'aitchison'}, optional
            Function to use when calculating the distance between two samples.
            Note that 'cityblock' and 'manhattan' are equivalent metrics.
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
        label : `string` or `callable`, optional
            A metadata field (or function) used to label each analysis. If passing a function, a
            dict containing the metadata for each analysis is passed as the first and only
            positional argument. The callable function must return a string.

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
        import altair as alt
        import numpy as np
        import pandas as pd
        from scipy.spatial.distance import squareform
        from scipy.stats import pearsonr
        from skbio.stats import ordination
        from sklearn import manifold
        from sklearn.metrics.pairwise import euclidean_distances

        if len(self._results) < 3:
            raise PlottingException(
                "There are too few samples for MDS/PCoA after filtering. Please select 3 or more "
                "samples to plot."
            )

        dists = self._compute_distance(rank, metric).to_data_frame()

        # here we figure out what to put in the tooltips and get the appropriate data
        if tooltip:
            if not isinstance(tooltip, list):
                tooltip = [tooltip]
        else:
            tooltip = []

        tooltip.insert(0, "Label")

        if color and color not in tooltip:
            tooltip.insert(1, color)

        if size and size not in tooltip:
            tooltip.insert(2, size)

        magic_metadata, magic_fields = self._metadata_fetch(tooltip, label=label)

        if method == OrdinationMethod.Smacof:
            # adapted from https://scikit-learn.org/stable/auto_examples/manifold/plot_mds.html
            x_field = "MDS1"
            y_field = "MDS2"

            seed = np.random.RandomState(seed=3)
            mds = manifold.MDS(
                max_iter=3000, eps=1e-12, random_state=seed, dissimilarity="precomputed", n_jobs=1
            )
            pos = mds.fit(dists).embedding_
            plot_data = pd.DataFrame(pos, columns=[x_field, y_field], index=dists.index)
            plot_data.index.names = ["classification_id"]
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
            x_extra_label = "r² = %.02f" % (r_squared[0],)
            y_extra_label = "r² = %.02f" % (r_squared[1],)
        elif method == OrdinationMethod.Pcoa:
            # suppress eigenvalue warning from skbio--not because it's an invalid warning, but
            # because lots of folks in the field run pcoa on these distances functions, even if
            # statistically inappropriate. perhaps this will change if we ever become more
            # opinionated about the analyses that we allow our users to do (roo)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ord_result = ordination.pcoa(
                    dists.round(6)
                )  # round to avoid float precision errors

            plot_data = ord_result.samples.iloc[:, [0, 1]]  # get first two components
            plot_data = plot_data.div(plot_data.abs().max(axis=0), axis=1)  # normalize to [0,1]
            plot_data.index = dists.index
            plot_data.index.names = ["classification_id"]
            x_field, y_field = plot_data.columns.tolist()  # name of first two components

            x_extra_label = "%0.02f%%" % (ord_result.proportion_explained[0] * 100,)
            y_extra_label = "%0.02f%%" % (ord_result.proportion_explained[1] * 100,)
        else:
            raise OneCodexException(
                "MDS method must be one of: {}".format(", ".join(OrdinationMethod.values))
            )

        # label the axes
        if xlabel is None:
            xlabel = "{} ({})".format(x_field, x_extra_label)
        if ylabel is None:
            ylabel = "{} ({})".format(y_field, y_extra_label)

        plot_data = pd.concat([plot_data, magic_metadata], axis=1).reset_index()
        plot_data["url"] = plot_data["classification_id"].apply(get_classification_url)

        alt_kwargs = dict(
            x=alt.X(x_field, axis=alt.Axis(title=xlabel)),
            y=alt.Y(y_field, axis=alt.Axis(title=ylabel)),
            tooltip=[magic_fields[t] for t in tooltip],
            href="url:N",
        )

        # only add these parameters if they are in use
        if color:
            color_kwargs = {
                "legend": alt.Legend(title=magic_fields[color]),
            }
            if not is_continuous(plot_data[color]) or has_missing_values(plot_data[color]):
                plot_data[color] = plot_data[color].fillna("N/A").astype(str)
                domain = plot_data[color].values
                color_range = interleave_palette(domain)
                color_kwargs["scale"] = alt.Scale(domain=domain, range=color_range)

            alt_kwargs["color"] = alt.Color(magic_fields[color], **color_kwargs)
        if size:
            alt_kwargs["size"] = magic_fields[size]

        chart = alt.Chart(plot_data).mark_circle(size=mark_size).encode(**alt_kwargs)

        chart = chart.properties(**prepare_props(title=title, height=height, width=width))
        open_links_in_new_tab(chart)

        if return_chart:
            return chart
        else:
            chart.interactive().display()
