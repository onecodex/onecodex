import warnings


class TaxonomyMixin(object):
    def tree_build(self):
        """Build a tree from the taxonomy data present in this `ClassificationsDataFrame` or
        `SampleCollection`.

        Returns
        -------
        `skbio.tree.TreeNode`, the root node of a tree that contains all the taxa in the current
        analysis and their parents leading back to the root node.
        """
        from skbio.tree import TreeNode

        # build all the nodes
        nodes = {}

        for tax_id in self.taxonomy.index:
            node = TreeNode(name=tax_id, length=1)
            node.tax_name = self.taxonomy["name"][tax_id]
            node.rank = self.taxonomy["rank"][tax_id]
            node.parent_tax_id = self.taxonomy["parent_tax_id"][tax_id]

            nodes[tax_id] = node

        # generate all the links
        for tax_id in self.taxonomy.index:
            try:
                parent = nodes[nodes[tax_id].parent_tax_id]
            except KeyError:
                if tax_id != "1":
                    warnings.warn(
                        "tax_id={} has parent_tax_id={} which is not in tree"
                        "".format(tax_id, nodes[tax_id].parent_tax_id)
                    )

                continue

            parent.append(nodes[tax_id])

        return nodes["1"]

    def tree_prune_tax_ids(self, tree, tax_ids):
        """Prunes a tree back to contain only the tax_ids in the list and their parents.

        Parameters
        ----------
        tree : `skbio.tree.TreeNode`
            The root node of the tree to perform this operation on.
        tax_ids : `list`
            A `list` of taxonomic IDs to keep in the tree.

        Returns
        -------
        `skbio.tree.TreeNode`, the root of a tree containing the given taxonomic IDs and their
        parents, leading back to the root node.
        """
        tax_ids_to_keep = []

        for tax_id in tax_ids:
            tax_ids_to_keep.append(tax_id)
            tax_ids_to_keep.extend([x.name for x in tree.find(tax_id).ancestors()])

        tree = tree.copy()
        tree.remove_deleted(lambda n: n.name not in tax_ids_to_keep)

        return tree

    def tree_prune_rank(self, tree, rank="species"):
        """Takes a TreeNode tree and prunes off any tips not at the specified rank and backwards up
        until all of the tips are at the specified rank.

        Parameters
        ----------
        tree : `skbio.tree.TreeNode`
            The root node of the tree to perform this operation on.
        rank : {kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.

        Returns
        -------
        `skbio.tree.TreeNode`, the root of the tree where all tips are at the given rank, and all
        tips have a path back to the root node.

        Examples
        --------
        tree_prune_rank(tree, 'species') will remove all subspecies/strain nodes and return a tree
        containing all genus-level nodes and higher.
        """
        if rank is None:
            return tree.copy()

        tree = tree.copy()

        for node in tree.postorder():
            if node.rank == rank:
                node._above_rank = True
            elif any([getattr(n, "_above_rank", False) for n in node.children]):
                node._above_rank = True
            else:
                node._above_rank = False

        tree.remove_deleted(lambda n: not getattr(n, "_above_rank", False))

        return tree
