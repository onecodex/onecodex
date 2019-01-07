from hashlib import sha256
import pandas as pd
import pytest
import warnings

import onecodex
from onecodex.exceptions import OneCodexException
from onecodex.helpers import ResultsDataFrame, ResultsSeries


# test ResultsDataFrame
# test ResultsSeries -is metadata transferred between DataFrame, Series, and back?
# a few tests for AnalysisMethods--magic_metadata_fetch, results filtering, guessing normalization
# a few tests for EnhancedSampleCollection--magic_classification_fetch, collate_metadata, collate_results

def test_pandas_subclass():
    inner_df = pd.DataFrame({'datum1': [7, 4, 21], 'datum2': [8, 16, 24]})

    ocx_data = {
        'ocx_rank': 'jedi',
        'ocx_field': 'bantha poodoo',
        'ocx_taxonomy': inner_df.copy(),
        'ocx_metadata': inner_df.copy()
    }

    df = ResultsDataFrame({'1279': [1, 2, 3], '1280': [4, 5, 6]}, ocx_data=ocx_data)

    # we want to be sure that slices of our df are returned as ResultsSeries
    assert type(df['1279']) is ResultsSeries

    # we're mostly interested in whether our metadata is transferred between copies. some operations
    # split the df into series and concat back to df, so by doing all this stuff to it we're actually
    # testing several consecutive copy operations using different parts of the pandas API
    new_df = (((df * 10) / 2.2).round(5)) ** 0.5

    # rank is explicitly /not/ passed on, since we don't know what the user has done to the df and
    # we therefore can't trust the rank to be correct
    assert new_df.ocx_rank is None
    assert new_df.ocx_field == 'bantha poodoo'
    assert (new_df.ocx_taxonomy == inner_df).all().all()
    assert (new_df.ocx_metadata == inner_df).all().all()


def test_auto_rank(ocx_w_enhanced, api_data):
    samples = ocx_w_enhanced.Samples.where(project='4b53797444f846c4')

    # if we re-collate the results using the abundance field, auto rank should choose species
    samples._collate_results(field='abundance')
    assert samples._get_auto_rank('auto') == 'species'

    # inside the pandas extension, auto rank should choose that of the ResultsDataFrame
    results = samples.results(rank='phylum')
    assert results.ocx._get_auto_rank('auto') == 'phylum'


def test_guess_normalization(ocx_w_enhanced, api_data):
    samples = ocx_w_enhanced.Samples.where(project='4b53797444f846c4')

    norm_results = samples.results(normalize=True)
    assert norm_results.ocx._guess_normalized() is True

    unnorm_results = samples.results(normalize=False)
    assert unnorm_results.ocx._guess_normalized() is False


def test_metadata_fetch(ocx_w_enhanced, api_data):
    samples = ocx_w_enhanced.Samples.where(project='4b53797444f846c4')

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


def test_results_filtering_rank(ocx_w_enhanced, api_data):
    samples = ocx_w_enhanced.Samples.where(project='4b53797444f846c4')

    tree = samples.tree_build()
    tree = samples.tree_prune_tax_ids(tree, ['1279'])

    # test filtering by rank
    genus_results = samples.results(rank='genus')
    family_results = samples.results(rank='family')

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
        samples.results(rank='does_not_exist')
    assert 'No taxa kept' in str(e.value)


def test_results_filtering_other(ocx_w_enhanced, api_data):
    samples = ocx_w_enhanced.Samples.where(project='4b53797444f846c4')

    # normalize the data
    norm_results = samples.results(normalize=True)
    assert norm_results.sum(axis=1).round(6).tolist() == [1.0, 1.0, 1.0]

    # trying to 'unnormalize' data fails
    with pytest.raises(OneCodexException) as e:
        norm_results.ocx.results(normalize=False)
    assert 'already been normalized' in str(e.value)

    # remove columns where every value is zero
    results = samples.results(rank=None, normalize=False, remove_zeros=False)
    assert len(results.columns) == 3157
    results['1279'] = 0
    results['1280'] = 0
    results = results.ocx.results(rank=None, normalize=False, remove_zeros=True)
    assert len(results.columns) == 3155

    # return only taxa with at least 100 reads in one or more samples
    assert ((samples.results(rank=None, normalize=False, remove_zeros=False, threshold=100) >= 100).any(axis=0)).all()

    # top N most abundant taxa
    assert sha256(samples.results(top_n=10).round(6).to_json().encode()).hexdigest() == \
        '437fcf282b885440571f552bf322056b6633154b247a9382ca074eee4b1ebf59'

    # check entire contents of long and wide format tables
    wide_tbl = samples.results(table_format='wide', normalize=False)
    assert wide_tbl.sum().sum() == 38034443
    assert sum(wide_tbl.columns.astype(int).tolist()) == 109170075

    long_tbl = samples.results(table_format='long', normalize=False)
    assert long_tbl['readcount_w_children'].sum() == 38034443
    assert long_tbl['tax_id'].astype(int).sum() / 3 == 109170075

    # should fail if format is not wide or long
    with pytest.raises(OneCodexException) as e:
        samples.results(table_format='does_not_exist')
    assert 'must be one of' in str(e.value)


def test_enhanced_sample_collection(ocx_w_enhanced, api_data):
    samples = ocx_w_enhanced.Samples.where(project='4b53797444f846c4')

    # manipulations of samples in the collection ought to update the stored dfs
    class_id = samples[2].primary_classification.id
    del samples[2]

    assert len(samples) == 2
    assert len(samples.results()) == 2
    assert len(samples.metadata) == 2
    assert class_id not in samples._results.index
    assert class_id not in samples.metadata.index


def test_classification_fetch(ocx_w_enhanced, api_data):
    samples = ocx_w_enhanced.Samples.where(project='4b53797444f846c4')

    # should work with a list of classifications as input, not just samples
    samples._oc_model = ocx_w_enhanced.Classifications
    samples._res_list = samples.primary_classifications
    samples._resource = [x._resource for x in samples.primary_classifications]
    samples._classification_fetch()

    # should issue a warning if a classification did not succeed
    with warnings.catch_warnings(record=True) as w:
        samples._resource[0].success = False
        samples._update()
        samples._classification_fetch()
        assert len(w) == 1
        assert 'not successful' in str(w[-1].message)

    # be sure to change success back to True, or other tests will ignore this classification
    # result--the resource object persists for the session!
    samples._resource[0].success = True


def test_collate_metadata(ocx_w_enhanced, api_data):
    samples = ocx_w_enhanced.Samples.where(project='4b53797444f846c4')

    # check contents of metadata df--at least that which can easily be coerced to strings
    metadata = samples.metadata
    string_to_hash = ""
    for col in sorted(metadata.columns.tolist()):
        for row in sorted(metadata.index.tolist()):
            try:
                string_to_hash += metadata.loc[row, col].astype(str)
            except AttributeError:
                pass

    assert sha256(string_to_hash.encode()).hexdigest() == \
        '6bebbdcc842f5d83d98a02657231093d68a649fc7721cf3a92755260dd45bf3d'

    # label must be a string or callable
    with pytest.raises(NotImplementedError) as e:
        samples._collate_metadata(label=123)
    assert 'string or function' in str(e.value)

    # label is a field not in the table
    with pytest.raises(OneCodexException) as e:
        samples._collate_metadata(label='does_not_exist')
    assert 'not find any labels' in str(e.value)


def test_collate_results(ocx_w_enhanced, api_data):
    samples = ocx_w_enhanced.Samples.where(project='4b53797444f846c4')

    # check contents of results df
    string_to_hash = ""
    for col in sorted(samples._results.columns.tolist()):
        for row in sorted(samples._results.index.tolist()):
            try:
                string_to_hash += samples._results.loc[row, col].astype(str)
            except AttributeError:
                pass

    assert sha256(string_to_hash.encode()).hexdigest() == \
        '5353e7ada1db22624db3604f80e3d6b2c15ca476c9d47881d0e8123b7e6c009d'

    # check contents of taxonomy df
    string_to_hash = ""
    for col in sorted(samples.taxonomy.columns.tolist()):
        for row in sorted(samples.taxonomy.index.tolist()):
            try:
                string_to_hash += samples.taxonomy.loc[row, col].astype(str)
            except AttributeError:
                pass

    assert sha256(string_to_hash.encode()).hexdigest() == \
        'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855'

    # invalid field name
    with pytest.raises(OneCodexException) as e:
        samples._collate_results(field='does_not_exist')
    assert 'not valid' in str(e.value)


def test_pandas_extension(ocx_w_enhanced, api_data):
    samples = ocx_w_enhanced.Samples.where(project='4b53797444f846c4')
    results = samples.results()

    # extension should be in ocx namespace of ResultsDataFrame
    assert getattr(results, 'ocx', False)
    assert isinstance(results.ocx, onecodex.helpers.OneCodexAccessor)
    assert type(results.ocx).__base__ == onecodex.helpers.AnalysisMethods

    # changes to contents of results df should affect contents of taxonomy df, by keeping only
    # tax_ids in the results df and their parents
    results = samples.results(top_n=2)
    assert sorted(results.ocx.taxonomy.index.tolist(), key=int) == \
        ['1', '2', '191', '815', '816', '976', '1224', '28211', '41295',
         '68336', '131567', '171549', '200643', '204441', '1783270']
