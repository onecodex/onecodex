import warnings


class TaxonomyMixin():
    def skbio_tree_build(self):
        from skbio.tree import TreeNode

        # build all the nodes
        nodes = {}

        for tax_id in self._taxonomy.index:
            node = TreeNode(name=tax_id, length=1)
            node.tax_name = self._taxonomy['name'][tax_id]
            node.rank = self._taxonomy['rank'][tax_id]
            node.parent_tax_id = self._taxonomy['parent_tax_id'][tax_id]

            nodes[tax_id] = node

        # generate all the links
        for tax_id in self._taxonomy.index:
            try:
                parent = nodes[nodes[tax_id].parent_tax_id]
            except KeyError:
                if tax_id != '1':
                    warnings.warn('tax_id={} has parent_tax_id={} which is not in tree'
                                  ''.format(tax_id, nodes[tax_id].parent_tax_id))

                continue

            parent.append(nodes[tax_id])

        return nodes['1']

    def skbio_tree_prune(self, tree, rank='species'):
        """
        Takes a TreeNode tree and prunes off any tips not at the specified rank
        and backwards up until all of the tips are at the specified rank.
        """

        tree = tree.copy()

        for node in tree.postorder():
            if node.rank == rank:
                node._above_rank = True
            elif any([getattr(n, '_above_rank', False) for n in node.children]):
                # TODO: should this ever happen?
                node._above_rank = True
            else:
                node._above_rank = False

        tree.remove_deleted(lambda n: not getattr(n, '_above_rank', False))

        return tree
