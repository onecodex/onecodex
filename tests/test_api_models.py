from __future__ import print_function
import datetime
import pytest
import responses
import sys

try:
    from urllib.parse import unquote_plus  # Py3
except ImportError:
    from urllib import unquote_plus

import onecodex
from onecodex import Api
from onecodex.exceptions import MethodNotSupported, OneCodexException


def test_api_creation(api_data):
    ocx = Api(api_key='1eab4217d30d42849dbde0cd1bb94e39',
              base_url='http://localhost:3000', cache_schema=False)
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


def test_sample_download(ocx, api_data):
    sample = ocx.Samples.get('761bc54b97f64980')

    sample.download()


def test_resourcelist(ocx, api_data):
    sample = ocx.Samples.get('761bc54b97f64980')
    tags1 = onecodex.models.ResourceList(sample.tags._resource,
                                         onecodex.models.misc.Tags)

    assert isinstance(sample.tags, onecodex.models.ResourceList)
    assert sample.tags == tags1

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
    assert id(tags1) != id(sample.tags)
    assert id(tags1._resource) == id(sample.tags._resource)

    # TODO: these tests shouldn't fail, but we're leaving this bug here for now due to conversations
    # with @boydgreenfield on 1/10/2019
    # assert len(tags1) == len(sample.tags)

    # for i in range(len(tags1)):
    #     assert tags1[i] == sample.tags[i]

    # can't mix types in a ResourceList
    with pytest.raises(ValueError) as e:
        tags1.append(sample)
    assert 'object of type' in str(e.value)

    if sys.version_info.major < 3:
        with pytest.raises(AttributeError):
            tags1.clear()
    else:
        tags1.clear()
        assert len(tags1) == 0


def test_samplecollection(ocx, api_data):
    all_samples = ocx.Samples.where()
    samples = all_samples[:3]
    other_samples = all_samples[4:7]

    # duplicate Classifications can not be part of the same SampleCollection
    with pytest.raises(OneCodexException) as e:
        samples + samples
    assert 'contain duplicate objects' in str(e.value)

    # SampleCollections can be added together
    new_samples = samples + other_samples
    assert len(samples) == len(other_samples) == 3
    assert len(new_samples) == 6

    # addition doesn't work with a SampleCollection and a lone Samples object
    single_sample = ocx.Samples.get('761bc54b97f64980')

    with pytest.raises(TypeError) as e:
        samples + single_sample
    assert 'can only concatenate' in str(e.value)


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

    counts = []

    for c in responses.calls:
        url = unquote_plus(c.request.url)

        # there may be accessory requests, like looking up a classification result associated
        # with a sample. therefore, we expect to see the query we're looking in for in one (but not
        # all) of the requests. the order of the requests is not deterministic, so check them all.
        count = 0

        for query in queries:
            if query in url:
                count += 1

        counts.append(count)

    # the correct queries must both appear together in only one request to pass this test
    assert len([x for x in counts if x == len(queries)]) == 1


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


def test_jobs(ocx, api_data):
    jobs = ocx.Jobs.all()
    assert len(jobs) == 23

    jobs = ocx.Jobs.where(public=True)
    assert len(jobs) == 23
