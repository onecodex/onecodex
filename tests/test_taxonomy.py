import pytest

from onecodex.lib.taxonomy import generate_skbio_tree, prune_to_rank


def test_tree_generation(ocx, api_data):
    classification = ocx.Classifications.get('f9e4a5506b154953')
    tree = generate_skbio_tree(classification)
    assert tree.has_children()
    assert tree.is_root()
    staph = tree.find('1279')
    assert staph.tax_name == 'Staphylococcus'


def test_tree_pruning(ocx, api_data):
    from skbio.tree import MissingNodeError
    classification = ocx.Classifications.get('f9e4a5506b154953')
    tree = generate_skbio_tree(classification)

    # the species and genus nodes exist
    tree.find('1279')
    tree.find('1078083')

    # prune back to genus
    prune_to_rank(tree, 'genus')

    # the species node no longer exists, but the genus still does
    with pytest.raises(MissingNodeError):
        tree.find('1078083')
    tree.find('1279')
