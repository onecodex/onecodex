from collections import defaultdict


def _merge_unlinked(tax_id, unlinked):
    # we might have to scan unlinked to see if any of the children have
    # unlinked nodes too (and so on recursively)
    children = unlinked.pop(tax_id, [])
    for child in children:
        child.extend(_merge_unlinked(child.name, unlinked))
    return children


def generate_skbio_tree(classification, existing_tree=None):
    from skbio.tree import MissingNodeError, TreeNode
    # we want to use NCBI tax ids and skbio randomly overwrites these
    TreeNode.assign_ids = lambda _: None

    otus = classification.results()['table']
    if existing_tree is None:
        tree = TreeNode(name='1')
        tree.tax_name = 'Root'
        tree.rank = 'no rank'
    else:
        tree = existing_tree

    # we use this to keep track of nodes that haven't had their parent added yet
    unlinked = defaultdict(list)

    for otu in otus:
        tax_id = otu['tax_id']
        # skip nodes already in the tree
        try:
            tree.find(tax_id)
            continue
        except MissingNodeError:
            pass

        # try to find a parent (if it exists)
        parent_id = otu['parent_tax_id']
        try:
            parent = tree.find(parent_id)
            # the children are merged out here (only if we have a parent) to
            # make sure we're not creating trees inside unlinked itself
            children = _merge_unlinked(tax_id, unlinked)
        except MissingNodeError:
            parent = None
            children = None

        # create the node
        node = TreeNode(name=tax_id, children=children)
        node.tax_name = otu.get('name', '')
        node.rank = otu.get('rank', 'no rank')

        # either add the node to its parent or keep track of it until its
        # parent is "in tree" too
        if parent is not None:
            parent.append(node)
        else:
            unlinked[parent_id].append(node)

    assert len(unlinked) == 0, 'some unlinked nodes were not included in the tree'

    return tree
