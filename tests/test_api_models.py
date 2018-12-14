from __future__ import print_function
import datetime
import json
import pandas as pd
import pytest
import responses

try:
    from urllib.parse import unquote_plus  # Py3
except ImportError:
    from urllib import unquote_plus

import onecodex
from onecodex import Api
from onecodex.exceptions import MethodNotSupported


def test_api_creation(api_data):
    ocx = Api(api_key='1eab4217d30d42849dbde0cd1bb94e39',
              base_url='http://localhost:3005', cache_schema=False)
    assert isinstance(ocx, Api)
    assert True


def test_sample_get(ocx, api_data):
    sample = ocx.Samples.get('761bc54b97f64980')
    assert sample.size == 302369471
    assert sample.filename == 'SRR2352223.fastq.gz'
    assert sample.__repr__() == '<Samples 761bc54b97f64980: "SRR2352223.fastq.gz">'
    assert isinstance(sample.created_at, datetime.datetime)

    classification = sample.primary_classification
    assert classification
    assert classification.complete

    tags = sample.tags
    assert len(tags) > 1
    assert 'isolate' in [t.name for t in tags]


def test_resourcelist(ocx, api_data):
    sample = ocx.Samples.get('761bc54b97f64980')
    tags1 = onecodex.models.ResourceList(sample.tags._resource,
                                         onecodex.models.misc.Tags)
    tags2 = onecodex.models.ResourceList(sample.tags._resource,
                                         onecodex.models.misc.Tags)

    assert isinstance(sample.tags, onecodex.models.ResourceList)

    # test manipulation of tags lists
    tag_to_pop = sample.tags[-1]
    popped_tag = sample.tags.pop()
    assert id(tag_to_pop._resource) == id(popped_tag._resource)

    sample.tags.insert(0, popped_tag)
    assert id(sample.tags[0]._resource) == id(popped_tag._resource)
    assert sample.tags.index(popped_tag) == 0

    assert sample.tags.count(popped_tag) == 1
    sample.tags.remove(popped_tag)
    assert sample.tags.count(popped_tag) == 0

    with pytest.raises(ValueError):
        sample.tags.remove(popped_tag)
    with pytest.raises(ValueError):
        sample.tags.index(popped_tag)

    # we can set tags list in-place
    sample.tags[0] = popped_tag
    assert id(sample.tags[0]._resource) == id(popped_tag._resource)

    # changes in one instance of a ResourceList affect other instances
    assert id(tags1) != id(tags2) != id(sample.tags)
    assert id(tags1._resource) == id(tags2._resource) == id(sample.tags._resource)

    assert len(tags1) == len(tags2) == len(sample.tags)

    for i in range(len(tags1)):
        assert tags1[i] == tags2[i] == sample.tags[i]


def test_comparisons(ocx, api_data):
    sample1 = ocx.Samples.get('761bc54b97f64980')
    sample2 = ocx.Samples.get('761bc54b97f64980')

    # should be different objects
    assert id(sample1) != id(sample2)

    # but they should reference the same Resource
    assert id(sample1._resource) == id(sample2._resource)

    # and so they should be equal
    assert sample1 == sample2

    # but any non-equivalence comparisons should throw errors
    with pytest.raises(NotImplementedError):
        sample1 > sample2
    with pytest.raises(NotImplementedError):
        sample1 >= sample2
    with pytest.raises(NotImplementedError):
        sample1 < sample2
    with pytest.raises(NotImplementedError):
        sample1 <= sample2

    # and this should be true of ResourceLists too
    assert sample1.tags == sample2.tags

    with pytest.raises(NotImplementedError):
        sample1.tags > sample2.tags
    with pytest.raises(NotImplementedError):
        sample1.tags >= sample2.tags
    with pytest.raises(NotImplementedError):
        sample1.tags < sample2.tags
    with pytest.raises(NotImplementedError):
        sample1.tags <= sample2.tags

    # and so should comparisons with asimilar objects
    assert sample1 != [1, 2, 3]
    assert sample1.tags != [1, 2, 3]


def test_dir_method(ocx, api_data):
    sample = ocx.Samples.get('761bc54b97f64980')

    instance_names = dir(sample)
    assert "where" not in instance_names   # we mask @classmethod's
    assert "created_at" in instance_names  # property on _resource
    assert "save" in instance_names        # function in py3, method in py2

    class_names = dir(ocx.Samples)
    assert "where" in class_names
    assert "created_at" not in class_names
    assert "save" in class_names  # instance methods are available off class


def test_get_failure_instructions(ocx):
    with pytest.raises(TypeError):
        ocx.Samples('direct_id')


def test_model_deletions(ocx, api_data):
    sample = ocx.Samples.get('761bc54b97f64980')
    sample.delete()

    analysis = sample.primary_classification
    with pytest.raises(MethodNotSupported):
        analysis.delete()


def test_model_updates(ocx, api_data):
    sample = ocx.Samples.get('761bc54b97f64980')
    sample.visibility = 'public' if sample.visibility == 'private' else 'public'

    # Read-only field
    with pytest.raises(MethodNotSupported):
        sample.filename = 'something_else'

    # No update resource
    analysis = sample.primary_classification
    with pytest.raises(MethodNotSupported):
        analysis.created_at = datetime.datetime.utcnow()


def test_sample_saving(ocx, api_data):
    sample = ocx.Samples.get('761bc54b97f64980')
    visibility = sample.visibility
    sample.visibility = 'private' if visibility == 'public' else 'public'
    sample.save()
    assert sample.visibility is not visibility


def test_metadata_saving(ocx, api_data):
    sample = ocx.Samples.get('761bc54b97f64980')
    metadata1 = sample.metadata
    metadata2 = ocx.Metadata.get('4fe05e748b5a4f0e')
    assert metadata1 == metadata2
    metadata1.description = 'my new description -- testing!'
    metadata1.date_collected = datetime.datetime.now()
    metadata1.save()
    assert isinstance(metadata1.date_collected, datetime.datetime)
    assert metadata1.description == 'my new description -- testing!'
    assert hasattr(metadata1, 'sample')  # This will fail because we don't mock it in the return


def test_dir_patching(ocx, api_data):
    sample = ocx.Samples.get('761bc54b97f64980')
    props = {'id', 'created_at', 'filename', 'visibility',
             'metadata', 'owner', 'primary_classification', 'project',
             'size', 'tags'}
    for prop in props:
        assert prop in dir(sample)
    assert len(sample.__dict__) == 1  # I'm not sure we *want* this...


def test_classification_methods(ocx, api_data):
    classification = ocx.Classifications.get('45a573fb7833449a')
    assert isinstance(classification, onecodex.models.analysis.Classifications)
    t = classification.table()
    assert isinstance(t, pd.DataFrame)


def test_no_results_on_generic_analysis(ocx, api_data):
    analysis = ocx.Analyses.get('45a573fb7833449a')
    with pytest.raises(NotImplementedError):
        analysis.results()


# Sorting and where clauses
@pytest.mark.parametrize('where_args,where_kwargs,queries', [
    ([], {'visibility': 'public'},
        ['where={"visibility": "public"}']),
    ([], {'visibility': 'private'},
        ['where={"visibility": "private"}']),
    ([], {'filename': 'SRR1234.fastq.gz'},
        ['where={"filename": "SRR1234.fastq.gz"}']),
    ([], {'filename': 'SRR1234.fastq.gz', 'sort': 'visibility'},
        ['where={"filename": "SRR1234.fastq.gz"}']),
    ([], {'visibility': 'private', 'filename': 'tmp.fa'},
        ['"filename": "tmp.fa"', '"visibility": "private"']),
    ([{'visibility': 'private', 'filename': 'tmp.fa'}], {},
        ['"filename": "tmp.fa"', '"visibility": "private"']),
    (['761bc54b97f64980'], {},
        ['"$uri": {"$in": ["/api/v1/samples/761bc54b97f64980"']),
    (['/api/v1/samples/761bc54b97f64980'], {},
        ['"$uri": {"$in": ["/api/v1/samples/761bc54b97f64980"']),
])
def test_where_clauses(ocx, api_data, where_args, where_kwargs, queries):
    ocx.Samples.where(*where_args, **where_kwargs)
    urls = []
    for c in responses.calls:
        url = unquote_plus(c.request.url)
        urls.append(url)
        for query in queries:
            assert query in url


def test_public_search(ocx, api_data):
    with pytest.warns(DeprecationWarning):
        samples = ocx.Samples.search_public(filename='tmp.fa')
        assert len(samples) == 0
    samples = ocx.Samples.where(filename='tmp.fa', public=True)
    assert len(samples) == 0


def test_public_project(ocx, api_data):
    with pytest.warns(DeprecationWarning):
        projs = ocx.Projects.search_public(name='One Codex Project')
        assert len(projs) == 0
    projs = ocx.Projects.where(name='One Codex Project', public=True)
    assert len(projs) == 0


def test_where_clauses_with_tags(ocx, api_data):
    tag = ocx.Tags.get('5c1e9e41043e4435')
    sample = ocx.Samples.get('7428cca4a3a04a8e')
    samples = ocx.Samples.where(tags=[tag])
    assert sample in samples

    query = '{"tags": {"$containsall": [{"$ref": "/api/v1/tags/5c1e9e41043e4435"}]}}'
    query_in_urls = []
    for c in responses.calls:
        url = unquote_plus(c.request.url)
        query_in_urls.append(query in url)

    assert any(query_in_urls)


def test_where_primary_classification(ocx, api_data):
    analysis = ocx.Analyses.get('935c2a3611944e39')
    sample = ocx.Samples.get('7428cca4a3a04a8e')
    samples = ocx.Samples.where(primary_classification=analysis)
    assert sample in samples

    query = '{"primary_classification": {"$ref": "/api/v1/analyses/935c2a3611944e39"}'
    query_in_urls = []
    for c in responses.calls:
        url = unquote_plus(c.request.url)
        query_in_urls.append(query in url)

    assert any(query_in_urls)


def test_public_analyses(ocx, api_data):
    analyses = ocx.Analyses.where(public=True)
    assert len(analyses) == 1
    a = analyses[0]
    assert a.sample.filename == 'MSA-1000.16S.example.fastq.gz'
    assert a.job.name == 'One Codex Database'
    assert a.sample.visibility == 'public'


def test_biom(ocx, api_data):
    c1 = ocx.Classifications.get('45a573fb7833449a')
    c2 = ocx.Classifications.get('593601a797914cbf')
    biom = ocx.Classifications.to_otu([c1, c2])
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


def test_jobs(ocx, api_data):
    jobs = ocx.Jobs.all()
    assert len(jobs) == 23

    jobs = ocx.Jobs.where(public=True)
    assert len(jobs) == 23


def test_modulealias():
    # ensure that viz and distance modules get imported into Api() instances
    ocx = Api()

    assert isinstance(ocx.viz, onecodex.utils.ModuleAlias)
    assert isinstance(ocx.distance, onecodex.utils.ModuleAlias)
