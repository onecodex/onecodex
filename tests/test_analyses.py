from hashlib import sha256
import pytest; pytest.importorskip('pandas')  # noqa
import warnings

from onecodex.exceptions import OneCodexException


def test_auto_rank(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')

    # if we re-collate the results using the abundance field, auto rank should choose species
    samples._collate_results(field='abundance')
    assert samples._get_auto_rank('auto') == 'species'

    # inside the pandas extension, auto rank should choose that of the ResultsDataFrame
    results = samples.to_df(rank='phylum')
    assert results.ocx._get_auto_rank('auto') == 'phylum'


def test_guess_normalization(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')

    norm_results = samples.to_df(normalize=True)
    assert norm_results.ocx._guess_normalized() is True

    unnorm_results = samples.to_df(normalize=False)
    assert unnorm_results.ocx._guess_normalized() is False


def test_metadata_fetch(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')

    # tuple with invalid field
    with pytest.raises(OneCodexException) as e:
        samples._metadata_fetch([('does_not_exist', 'vegetables')])
    assert 'not found' in str(e.value)

    # tuple where one field is not categorical
    with pytest.raises(OneCodexException) as e:
        samples._metadata_fetch([('totalige', 'vegetables')])
    assert 'must be categorical' in str(e.value)

    # proper tuple of categorical metadata fields
    df, fields = samples._metadata_fetch([('eggs', 'vegetables')])
    assert fields[('eggs', 'vegetables')] == 'eggs_vegetables'
    assert df['eggs_vegetables'].tolist() == ['True_True', 'True_True', 'True_True']

    # Label reserved field name
    df, fields = samples._metadata_fetch(['Label'])
    assert fields['Label'] == 'Label'
    assert df['Label'].tolist() == ['SRR4408293.fastq', 'SRR4305031.fastq', 'SRR4408292.fastq']

    # exact match of single metadata field
    df, fields = samples._metadata_fetch(['totalige'])
    assert fields['totalige'] == 'totalige'
    assert df['totalige'].tolist() == [62.9, 91.5, 112.0]

    # tax_id coerced to string from integer
    df, fields = samples._metadata_fetch([1279])
    assert fields[1279] == 'Staphylococcus (1279)'
    assert df['Staphylococcus (1279)'].round(10).tolist() == [4.0172e-06, 4.89491e-05, 2.0881e-06]

    # tax_name, should take lowest matching tax_id, be case-insensitive, and report within-rank abundance
    df, fields = samples._metadata_fetch(['bacteroid'])
    assert fields['bacteroid'] == 'Bacteroidaceae (815)'
    assert df['Bacteroidaceae (815)'].round(10).tolist() == [0.344903898, 0.1656058794, 0.7776093433]


def test_results_filtering_rank(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')

    tree = samples.tree_build()
    tree = samples.tree_prune_tax_ids(tree, ['1279'])

    # test filtering by rank
    genus_results = samples.to_df(rank='genus')
    family_results = samples.to_df(rank='family')

    assert '1279' in genus_results
    assert '1279' not in family_results
    assert '90964' in family_results
    assert '90964' not in genus_results

    # the resulting taxonomy dicts should have taxa at rank specified and their parents,
    # but no children, and no higher rank nodes not connected to taxa at the specified rank
    assert len(genus_results.ocx.taxonomy.loc[genus_results.ocx.taxonomy['rank'] == 'species']) == 0
    assert len(genus_results.ocx.taxonomy.loc[genus_results.ocx.taxonomy['rank'] == 'genus']) == 474
    assert len(genus_results.ocx.taxonomy.loc[genus_results.ocx.taxonomy['rank'] == 'family']) == 163
    assert len(family_results.ocx.taxonomy.loc[family_results.ocx.taxonomy['rank'] == 'species']) == 0
    assert len(family_results.ocx.taxonomy.loc[family_results.ocx.taxonomy['rank'] == 'genus']) == 0
    assert len(family_results.ocx.taxonomy.loc[family_results.ocx.taxonomy['rank'] == 'family']) == 177

    # catch invalid ranks
    with pytest.raises(OneCodexException) as e:
        samples.to_df(rank='does_not_exist')
    assert 'No taxa kept' in str(e.value)

    # should warn if using rank=kingdom
    with warnings.catch_warnings(record=True) as w:
        samples.to_df(rank='kingdom')
        assert len(w) == 1
        assert 'superkingdom' in str(w[-1].message)


def test_results_filtering_other(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')

    # normalize the data
    norm_results = samples.to_df(normalize=True)
    assert norm_results.sum(axis=1).round(6).tolist() == [1.0, 1.0, 1.0]

    # trying to 'unnormalize' data fails
    with pytest.raises(OneCodexException) as e:
        norm_results.ocx.to_df(normalize=False)
    assert 'already been normalized' in str(e.value)

    # remove columns where every value is zero
    results = samples.to_df(rank=None, normalize=False, remove_zeros=False)
    assert len(results.columns) == 3157
    results['1279'] = 0
    results['1280'] = 0
    results = results.ocx.to_df(rank=None, normalize=False, remove_zeros=True)
    assert len(results.columns) == 3155

    # return only taxa with at least 100 reads in one or more samples
    assert ((samples.to_df(rank=None, normalize=False, remove_zeros=False, threshold=100) >= 100).any(axis=0)).all()

    # top N most abundant taxa
    assert sha256(samples.to_df(top_n=10).round(6).to_json().encode()).hexdigest() == \
        '437fcf282b885440571f552bf322056b6633154b247a9382ca074eee4b1ebf59'

    # check entire contents of long and wide format tables
    wide_tbl = samples.to_df(table_format='wide', normalize=False)
    assert wide_tbl.sum().sum() == 38034443
    assert sum(wide_tbl.columns.astype(int).tolist()) == 109170075

    long_tbl = samples.to_df(table_format='long', normalize=False)
    assert long_tbl['readcount_w_children'].sum() == 38034443
    assert long_tbl['tax_id'].astype(int).sum() / 3 == 109170075

    # should fail if format is not wide or long
    with pytest.raises(OneCodexException) as e:
        samples.to_df(table_format='does_not_exist')
    assert 'must be one of' in str(e.value)
