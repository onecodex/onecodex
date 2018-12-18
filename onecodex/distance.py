import pandas as pd
import skbio.diversity
import warnings

from onecodex.exceptions import OneCodexException
from onecodex.taxonomy import TaxonomyMixin


class DistanceMixin(TaxonomyMixin):
    def alpha_diversity(self, metric='simpson', rank='species'):
        """Caculate the diversity within a community"""

        if metric not in ('simpson', 'chao1'):
            raise OneCodexException('For alpha diversity, metric must be one of: simpson, chao1')

        df = self.results(rank=rank, normalize=False)

        output = {
            'classification_id': [],
            metric: []
        }

        for c_id in df.index:
            output['classification_id'].append(c_id)
            output[metric].append(skbio.diversity.alpha_diversity(
                metric,
                df.loc[c_id].tolist(),
                [c_id]
            ).values[0])

        return pd.DataFrame(output).set_index('classification_id')

    def beta_diversity(self, metric='braycurtis', rank='species'):
        """Calculate the diversity between two communities"""

        if metric not in ('jaccard', 'braycurtis', 'cityblock'):
            raise OneCodexException('For beta diversity, metric must be one of: jaccard, braycurtis, cityblock')

        # unifrac needs read counts, not relative abundances
        if self.field == 'abundance':
            raise OneCodexException('Beta diversity requires unnormalized read counts. Use field=readcount_w_children')

        df = self.results(rank=rank, normalize=False)

        counts = []
        for c_id in df.index:
            counts.append(df.loc[c_id].tolist())

        return skbio.diversity.beta_diversity(metric, counts, df.index.tolist())

    def unifrac(self, weighted=True, strict=False, rank='species'):
        """
        A beta diversity metric that takes into account the relative relatedness of community members.
        Weighted UniFrac looks at abundances, unweighted UniFrac looks at presence
        """

        # unifrac needs read counts, not relative abundances
        if self.field == 'abundance':
            raise OneCodexException('UniFrac requires unnormalized read counts. Use field=readcount_w_children')

        df = self.results(rank=rank, normalize=False)

        counts = []
        for c_id in df.index:
            counts.append(df.loc[c_id].tolist())

        tax_ids = df.keys().tolist()

        tree = self.skbio_tree_build()
        tree = self.skbio_tree_prune(tree, rank=rank)

        # there's a bug (?) in skbio where it expects the root to only have
        # one child, so we do a little faking here
        from skbio.tree import TreeNode
        new_tree = TreeNode(name='fake root')
        new_tree.rank = 'no rank'
        new_tree.append(tree)

        # find nodes that aren't tips, prune, and warn--they may be errors in taxonomy
        counts_to_remove = []

        for t_id in tax_ids:
            node = tree.find(t_id)

            if not node.is_tip():
                warnings.warn('Found tax_id={} as an internal node (children={}), expected a tip'
                              ''.format(t_id, ', '.join([x.name for x in node.children])))
            else:
                continue

            # remove its children, keep the original node
            for child in node.children:
                counts_to_remove.append(tax_ids.index(child.name))
                tax_ids.remove(child.name)

            node.children = []

        for count_idx in counts_to_remove:
            for sample_idx in range(len(counts)):
                del counts[sample_idx][count_idx]

        # then finally run the calculation and return
        if weighted:
            return skbio.diversity.beta_diversity('weighted_unifrac', counts, df.index.tolist(),
                                                  tree=new_tree, otu_ids=tax_ids)
        else:
            return skbio.diversity.beta_diversity('unweighted_unifrac', counts, df.index.tolist(),
                                                  tree=new_tree, otu_ids=tax_ids)
