from hashlib import sha256
import json
import pytest; pytest.importorskip('pandas')  # noqa
import warnings

from onecodex.models.collection import SampleCollection
from onecodex.exceptions import OneCodexException


def test_sample_collection_pandas(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')

    # manipulations of samples in the collection ought to update the stored dfs
    class_id = samples[2].primary_classification.id
    del samples[2]

    assert len(samples) == 2
    assert len(samples.to_df()) == 2
    assert len(samples.metadata) == 2
    assert class_id not in samples._results.index
    assert class_id not in samples.metadata.index


def test_biom(ocx, api_data):
    c1 = ocx.Classifications.get('45a573fb7833449a')._resource
    c2 = ocx.Classifications.get('593601a797914cbf')._resource
    biom = SampleCollection([c1, c2], ocx.Classifications).to_otu()
    assert set(biom.keys()) == {
        'columns', 'data', 'date', 'format', 'format_url',
        'generated_by', 'id', 'matrix_element_type', 'matrix_type', 'rows', 'shape', 'type'
    }
    assert biom['format_url'] == 'http://biom-format.org'
    assert set(biom['columns'][0].keys()) == {'id', 'sample_id', 'sample_filename', 'metadata'}
    assert set(biom['rows'][0].keys()) == {'id', 'metadata'}
    assert 'taxonomy' in biom['rows'][0]['metadata']

    # IDs
    assert biom['columns'][0]['id'] == c1.id
    assert biom['columns'][0]['sample_id'] == c1.sample.id

    # Reults
    assert set(row['metadata']['taxonomy'] for row in biom['rows']) == {
        'Staphylococcus sp. HGB0015', 'Staphylococcus'
    }

    # Format is row_id, sample id (column), count
    assert biom['data'][0] == [0, 0, 3]
    assert biom['data'][1] == [0, 1, 80]
    assert biom['data'][2] == [1, 0, 0]
    assert biom['data'][3] == [1, 1, 0]
    assert biom['rows'][0]['id'] == '1078083'
    assert biom['rows'][1]['id'] == '1279'

    # Check that we're not including unnecessary references
    assert '$uri' not in biom['columns'][0]['metadata']
    assert 'sample' not in biom['columns'][0]['metadata']

    # Test serialization
    assert json.loads(json.dumps(biom)) == biom   # tests json serialization


def test_classification_fetch(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')

    # should work with a list of classifications as input, not just samples
    samples._oc_model = ocx.Classifications
    samples._res_list = samples._classifications
    samples._resource = [x._resource for x in samples._classifications]
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


def test_collate_metadata(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')

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


def test_collate_results(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')

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
