import warnings

from onecodex.exceptions import OneCodexException
from onecodex.taxonomy import TaxonomyMixin
from onecodex.lib.enums import AlphaDiversityMetric, BetaDiversityMetric, Rank


class DistanceMixin(TaxonomyMixin):
    def alpha_diversity(self, metric=AlphaDiversityMetric.Shannon, rank=Rank.Auto):
        """Calculate the diversity within a community.

        Parameters
        ----------
        metric : {'simpson', 'observed_taxa', 'shannon'}
            The diversity metric to calculate.
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.

        Returns
        -------
        pandas.DataFrame, a distance matrix.
        """
        import pandas as pd
        import skbio.diversity

        if not AlphaDiversityMetric.has_value(metric):
            raise OneCodexException(
                "For alpha diversity, metric must be one of: {}".format(
                    ", ".join(AlphaDiversityMetric.values())
                )
            )

        if metric == "chao1":
            warnings.warn(
                "`Chao1` is deprecated and will be removed in a future release. Please use `observed_taxa` instead.",
                DeprecationWarning,
            )

        df = self.to_df(rank=rank, normalize=self._guess_normalized())

        skbio_metric = "observed_otus" if metric == "observed_taxa" else metric
        output = skbio.diversity.alpha_diversity(skbio_metric, df.values, df.index, validate=False)

        return pd.DataFrame(output, columns=[metric])

    def beta_diversity(self, metric=BetaDiversityMetric.BrayCurtis, rank=Rank.Auto):
        """Calculate the diversity between two communities.

        Parameters
        ----------
        metric : {'jaccard', 'braycurtis', 'cityblock'}
            The distance metric to calculate.
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.

        Returns
        -------
        skbio.stats.distance.DistanceMatrix, a distance matrix.
        """
        import skbio.diversity

        if not BetaDiversityMetric.has_value(metric):
            raise OneCodexException(
                "For beta diversity, metric must be one of: {}".format(
                    ", ".join(BetaDiversityMetric.values())
                )
            )

        df = self.to_df(rank=rank, normalize=self._guess_normalized())

        if metric == BetaDiversityMetric.WeightedUnifrac:
            return self.unifrac(weighted=True, rank=rank)
        elif metric == BetaDiversityMetric.UnweightedUnifrac:
            return self.unifrac(weighted=False, rank=rank)
        elif metric == BetaDiversityMetric.Jaccard:
            df = df > 0  # Jaccard requires a boolean matrix, otherwise it throws a warning

        # NOTE: see #291 for a discussion on using these metrics with normalized read counts. we are
        # explicitly disabling skbio's check for a counts matrix to allow normalized data to make
        # its way into this function.
        return skbio.diversity.beta_diversity(metric, df.values, df.index, validate=False)

    def unifrac(self, weighted=True, rank=Rank.Auto):
        """Calculate the UniFrac beta diversity metric.

        UniFrac takes into account the relatedness of community members. Weighted UniFrac considers
        abundances, unweighted UniFrac considers presence.

        Parameters
        ----------
        weighted : `bool`
            Calculate the weighted (True) or unweighted (False) distance metric.
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.

        Returns
        -------
        skbio.stats.distance.DistanceMatrix, a distance matrix.
        """
        import skbio.diversity

        df = self.to_df(rank=rank, normalize=self._guess_normalized())

        ocx_rank = df.ocx_rank
        # The scikit-bio implementations of phylogenetic metrics require integer counts
        if self._guess_normalized():
            df = df * 10e9

        tax_ids = df.keys().tolist()

        tree = self.tree_build()
        tree = self.tree_prune_rank(tree, rank=ocx_rank)

        # `scikit-bio` requires that the tree root has no more than 2
        # children, otherwise it considers it "unrooted".
        #
        # https://github.com/biocore/scikit-bio/blob/f3ae1dcfe8ea88e52e19f6693d79e529d05bda04/skbio/diversity/_util.py#L89
        #
        # Our taxonomy root regularly has more than 2 children, so we
        # add a fake parent of `root` to the tree here.
        from skbio.tree import TreeNode

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
                otu_ids=tax_ids,
                normalized=True,
            )
        else:
            return skbio.diversity.beta_diversity(
                BetaDiversityMetric.UnweightedUnifrac, df, df.index, tree=new_tree, otu_ids=tax_ids,
            )
