from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

from onecodex.exceptions import OneCodexException
from onecodex.lib.enums import (
    AlphaDiversityMetric,
    BetaDiversityMetric,
    Metric,
    Rank,
)
from onecodex.models.base_sample_collection import BaseSampleCollection
from onecodex.taxonomy import TaxonomyMixin

if TYPE_CHECKING:
    import pandas as pd
    import skbio.stats.distance.DistanceMatrix
    import sklearn.metrics.DistanceMetric


class DistanceMixin(BaseSampleCollection, TaxonomyMixin):
    def alpha_diversity(
        self,
        metric: Metric = Metric.Auto,
        rank=Rank.Auto,
        diversity_metric: AlphaDiversityMetric = AlphaDiversityMetric.Shannon,
    ) -> pd.DataFrame:
        """Calculate the diversity within a community.

        Parameters
        ----------
        metric: Metric
            The taxonomic abundance metric to use. See onecodex.lib.enums.Metric for definitions.
        diversity_metric : AlphaDiversityMetric
            The diversity metric to calculate. Note that Shannon diversity is calculated using log
            base ``e`` (natural log).
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.


        Returns
        -------
        pandas.DataFrame, a distance matrix.
        """
        import pandas as pd
        import skbio.diversity

        if not AlphaDiversityMetric.has_value(diversity_metric):
            raise OneCodexException(
                "For alpha diversity, metric must be one of: {}".format(
                    ", ".join(AlphaDiversityMetric.values())
                )
            )

        if diversity_metric == AlphaDiversityMetric.Chao1:
            warnings.warn(
                "`Chao1` is deprecated and will be removed in a future release. Please use `observed_taxa` instead.",
                DeprecationWarning,
            )

        metric, rank = self._parse_classification_config_args(metric=metric, rank=rank)
        df = self.to_df(rank=rank, metric=metric)

        # special-case mapping to skbio value
        skbio_metric = (
            "sobs" if diversity_metric == AlphaDiversityMetric.ObservedTaxa else diversity_metric
        )
        output = skbio.diversity.alpha_diversity(skbio_metric, df.values, df.index, validate=False)

        return pd.DataFrame(output, columns=[diversity_metric])

    def beta_diversity(
        self,
        metric: Metric = Metric.Auto,
        rank: Rank = Rank.Auto,
        diversity_metric: BetaDiversityMetric = BetaDiversityMetric.BrayCurtis,
    ) -> skbio.stats.distance.DistanceMatrix:
        """Calculate the diversity between two communities.

        Parameters
        ----------
        metric: Metric
            The taxonomic abundance metric to use. See onecodex.lib.enums.Metric for definitions.
        diversity_metric : BetaDiversityMetric
            The distance metric to calculate.
            Note that 'cityblock' and 'manhattan' are equivalent metrics.
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.

        Returns
        -------
        skbio.stats.distance.DistanceMatrix, a distance matrix.
        """
        import skbio.diversity

        if not BetaDiversityMetric.has_value(diversity_metric):
            raise OneCodexException(
                "For beta diversity, metric must be one of: {}".format(
                    ", ".join(BetaDiversityMetric.values())
                )
            )

        if diversity_metric == BetaDiversityMetric.WeightedUnifrac:
            return self.unifrac(weighted=True, metric=metric, rank=rank)
        elif diversity_metric == BetaDiversityMetric.UnweightedUnifrac:
            return self.unifrac(weighted=False, metric=metric, rank=rank)
        elif diversity_metric == BetaDiversityMetric.Aitchison:
            return self.aitchison_distance(metric=metric, rank=rank)

        metric, rank = self._parse_classification_config_args(metric=metric, rank=rank)
        df = self.to_df(rank=rank, metric=metric)

        if diversity_metric == BetaDiversityMetric.Jaccard:
            df = df > 0  # Jaccard requires a boolean matrix, otherwise it throws a warning

        # NOTE: see #291 for a discussion on using these metrics with normalized read counts. we are
        # explicitly disabling skbio's check for a counts matrix to allow normalized data to make
        # its way into this function.
        skbio_metric = (
            "cityblock" if diversity_metric == BetaDiversityMetric.Manhattan else diversity_metric
        )
        return skbio.diversity.beta_diversity(
            skbio_metric,
            df.values,
            df.index,
            validate=False,
            pairwise_func=self._pairwise_distances,
        )

    def _pairwise_distances(
        self,
        abundances,
        # must be called metric and not diversity_metric (this gets called by skbio)
        metric: sklearn.metrics.DistanceMetric,
    ) -> skbio.stats.distance.DistanceMatrix:
        import numpy as np
        from scipy.spatial.distance import pdist, squareform

        # We're using scipy's `pdist` to compute pairwise distances because it only computes one
        # triangle of distances. `pdist` used to be scikit-bio's default, but
        # `sklearn.metrics.pairwise_distances` is the default in more recent versions of scikit-bio.
        # sklearn's `pairwise_distances` function computes both lower and upper triangles and can
        # result in an asymmetric matrix due to floating-point/roundoff error, which causes
        # `skbio.DistanceMatrix` to raise an error about symmetry. We've seen this happen to users
        # when computing cityblock distances from normalized abundances.
        distance_matrix = squareform(
            pdist(abundances, metric=metric), force="tomatrix", checks=False
        )

        # Bray-Curtis can return nan distances when comparing samples with all zero abundances (e.g.
        # [0, 0, 0] vs [0, 0, 0]). We replace nan distances with zero because:
        #
        # 1) this is a reasonable substitute value (R packages do this by default)
        # 2) it matches how the other distance metrics behave in this case
        # 3) it prevents `skbio.DistanceMatrix` from erroring out on nan values
        nan_distances = np.isnan(distance_matrix)
        if (
            metric == BetaDiversityMetric.BrayCurtis
            and nan_distances.any()
            and (abundances >= 0.0).all()
        ):
            distance_matrix[nan_distances] = 0.0

        return distance_matrix

    def unifrac(self, metric: Metric = Metric.Auto, rank: Rank = Rank.Auto, weighted: bool = True):
        """Calculate the UniFrac beta diversity metric.

        UniFrac takes into account the relatedness of community members. Weighted UniFrac considers
        abundances, unweighted UniFrac considers presence.

        Parameters
        ----------
        metric: Metric
            the taxonomic abundance metric to use. See onecodex.lib.enums.Metric for definitions
        weighted : `bool`
            Calculate the weighted (True) or unweighted (False) distance metric.
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.

        Returns
        -------
        skbio.stats.distance.DistanceMatrix, a distance matrix.
        """
        import skbio.diversity
        from skbio.tree import TreeNode

        metric, rank = self._parse_classification_config_args(metric=metric, rank=rank)
        df = self.to_df(rank=rank, metric=metric, fill_missing=True)

        # The scikit-bio implementations of phylogenetic metrics require integer counts
        if metric.is_normalized:
            df = df * 10e9

        tax_ids = df.keys().tolist()

        tree = self.tree_build()
        tree = self.tree_prune_rank(tree, rank=rank)

        # `scikit-bio` requires that the tree root has no more than 2
        # children, otherwise it considers it "unrooted".
        #
        # https://github.com/biocore/scikit-bio/blob/f3ae1dcfe8ea88e52e19f6693d79e529d05bda04/skbio/diversity/_util.py#L89
        #
        # Our taxonomy root regularly has more than 2 children, so we
        # add a fake parent of `root` to the tree here.

        new_tree = TreeNode(name="fake root")
        new_tree.rank = "no rank"
        new_tree.append(tree)

        # then finally run the calculation and return
        if weighted:
            return skbio.diversity.beta_diversity(
                BetaDiversityMetric.WeightedUnifrac,
                df,
                df.index,
                tree=new_tree,
                taxa=tax_ids,
                normalized=True,
            )
        else:
            return skbio.diversity.beta_diversity(
                BetaDiversityMetric.UnweightedUnifrac,
                df,
                df.index,
                tree=new_tree,
                taxa=tax_ids,
            )

    def aitchison_distance(
        self, metric: Metric = Metric.Auto, rank: Rank = Rank.Auto
    ) -> skbio.stats.distance.DistanceMatrix:
        """Calculate the Aitchison distance between samples.

        Aitchison distance is the Euclidean distance between centre logratio-normalized samples (abundances).
        As this requires log-transforms, we first need to 'estimate' zeros in the data;
        i.e. replace zeros with small, positive values, while maintaining a constant sum to 1.

        Parameters
        ----------
        metric: Metric
            the taxonomic abundance metric to use. See onecodex.lib.enums.Metric for definitions
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.

        Returns
        -------
        skbio.stats.distance.DistanceMatrix, a distance matrix.
        """
        import numpy as np
        from skbio.stats.composition import clr, multi_replace
        from skbio.stats.distance import DistanceMatrix
        from sklearn.metrics.pairwise import euclidean_distances

        metric, rank = self._parse_classification_config_args(metric=metric, rank=rank)
        df = self.to_df(metric=metric, rank=rank)

        df_n0 = multi_replace(df)  # replace 0s with positive small numbers
        df_n0_clr = clr(df_n0)  # clr-normalize
        aitchison_array = euclidean_distances(df_n0_clr, df_n0_clr)  # get the euclidean distances

        # Due to rounding differences, we must force mirroring on the matrix
        aitchison_dm = np.zeros(aitchison_array.shape)
        aitchison_dm[np.triu_indices(aitchison_array.shape[0], k=0)] = aitchison_array[
            np.triu_indices(aitchison_array.shape[0], k=0)
        ]
        aitchison_dm = aitchison_dm + aitchison_dm.T - np.diag(np.diag(aitchison_dm))
        aitchison_dm = DistanceMatrix(aitchison_dm, df.index)

        return aitchison_dm
