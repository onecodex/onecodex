from onecodex.lib.taxonomy import generate_skbio_tree


def test_tree_generation(ocx, api_data):
    classification = ocx.Classifications.get('f9e4a5506b154953')
    tree = generate_skbio_tree(classification)
    assert tree.has_children()
    assert tree.is_root()
    staph = tree.find('1279')
    assert staph.tax_name == 'Staphylococcus'
