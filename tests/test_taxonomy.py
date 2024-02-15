import pytest

pytest.importorskip("pandas")  # noqa
from skbio.tree import MissingNodeError


def test_tree_generation(samples):
    tree = samples.tree_build()
    assert tree.has_children()
    assert tree.is_root()
    staph = tree.find("1279")
    assert staph.tax_name == "Staphylococcus"


def test_tree_pruning(samples):
    tree = samples.tree_build()

    # the species and genus nodes exist
    tree.find("1279")
    tree.find("1280")

    # prune back to genus
    pruned_tree = samples.tree_prune_rank(tree, "genus")

    # the species node no longer exists, but the genus still does
    with pytest.raises(MissingNodeError):
        pruned_tree.find("1280")
    pruned_tree.find("1279")

    # prune back to only 1279 and its parents
    pruned_tree = samples.tree_prune_tax_ids(tree, ["1279"])

    for tax_id in ["1", "131567", "2", "1783272", "1239", "91061", "1385", "90964", "1279"]:
        pruned_tree.find(tax_id)
