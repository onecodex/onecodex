from __future__ import print_function
import datetime
import pandas as pd

import onecodex
from onecodex import Api
from onecodex.exceptions import MethodNotSupported

import pytest
import responses

try:
    from urllib.parse import unquote_plus  # Py3
except ImportError:
    from urllib import unquote_plus


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

    analysis = sample.primary_analysis
    assert analysis
    assert analysis.complete

    tags = sample.tags
    assert len(tags) > 1
    assert 'isolate' in [t.name for t in tags]


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

    analysis = sample.primary_analysis
    with pytest.raises(MethodNotSupported):
        analysis.delete()


def test_model_updates(ocx, api_data):
    sample = ocx.Samples.get('761bc54b97f64980')
    sample.starred = not sample.starred

    # Read-only field
    with pytest.raises(MethodNotSupported):
        sample.public = not sample.public

    # No update resource
    analysis = sample.primary_analysis
    with pytest.raises(MethodNotSupported):
        analysis.created_at = datetime.datetime.utcnow()


def test_sample_saving(ocx, api_data):
    sample = ocx.Samples.get('761bc54b97f64980')
    starred_state = sample.starred
    sample.starred = not starred_state
    sample.save()
    assert sample.starred is not starred_state


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
    props = {'id', 'created_at', 'filename', 'indexed', 'public',
             'metadata', 'owner', 'primary_analysis', 'project',
             'size', 'starred', 'tags'}
    for prop in props:
        assert prop in dir(sample)
    assert len(sample.__dict__) == 1  # I'm not sure we *want* this...


def test_classification_methods(ocx, api_data):
    classification = ocx.Classifications.get('f9e4a5506b154953')
    assert isinstance(classification, onecodex.models.analysis.Classifications)
    t = classification.table()
    assert isinstance(t, pd.DataFrame)


def test_no_results_on_generic_analysis(ocx, api_data):
    analysis = ocx.Analyses.get('f9e4a5506b154953')
    with pytest.raises(NotImplementedError):
        analysis.results()


# Sorting and where clauses
@pytest.mark.parametrize('where_args,where_kwargs,queries', [
    ([], {'public': True},
        ['where={"public": true}']),
    ([], {'public': False},
        ['where={"public": false}']),
    ([], {'filename': 'SRR1234.fastq.gz'},
        ['where={"filename": "SRR1234.fastq.gz"}']),
    ([], {'filename': 'SRR1234.fastq.gz', 'sort': 'public'},
        ['where={"filename": "SRR1234.fastq.gz"}']),
    ([], {'public': False, 'filename': 'tmp.fa'},
        ['"filename": "tmp.fa"', '"public": false']),
    ([{'public': False, 'filename': 'tmp.fa'}], {},
        ['"filename": "tmp.fa"', '"public": false']),
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
    samples = ocx.Samples.search_public(filename="tmp.fa")
    assert len(samples) == 0


def test_public_project(ocx, api_data):
    projs = ocx.Projects.search_public(name="One Codex Project")
    assert len(projs) == 0


def test_where_clauses_with_tags(ocx, api_data):
    tag = ocx.Tags.get('5c1e9e41043e4435')
    sample = ocx.Samples.get('761bc54b97f64980')
    samples = ocx.Samples.where(tags=[tag])
    assert sample in samples

    query = '{"tags": {"$containsall": [{"$ref": "/api/v1/tags/5c1e9e41043e4435"}]}}'
    query_in_urls = []
    for c in responses.calls:
        url = unquote_plus(c.request.url)
        query_in_urls.append(query in url)

    assert any(query_in_urls)


def test_where_primary_analysis(ocx, api_data):
    analysis = ocx.Analyses.get('935c2a3611944e39')
    sample = ocx.Samples.get('761bc54b97f64980')
    samples = ocx.Samples.where(primary_analysis=analysis)
    assert sample in samples

    query = '{"primary_analysis": {"$ref": "/api/v1/analyses/935c2a3611944e39"}'
    query_in_urls = []
    for c in responses.calls:
        url = unquote_plus(c.request.url)
        query_in_urls.append(query in url)

    assert any(query_in_urls)
