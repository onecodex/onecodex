from __future__ import print_function
import datetime
import io
import json
import mock
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
    ocx = Api(
        api_key="1eab4217d30d42849dbde0cd1bb94e39",
        base_url="http://localhost:3000",
        cache_schema=False,
    )
    assert isinstance(ocx, Api)


def test_retries_set_on_client_session(api_data):
    ocx = Api(
        api_key="1eab4217d30d42849dbde0cd1bb94e39",
        base_url="http://localhost:3000",
        cache_schema=False,
    )

    assert ocx._session.adapters["http://"].max_retries.total == 3
    assert ocx._session.adapters["http://"].max_retries.allowed_methods is None
    assert ocx._session.adapters["https://"].max_retries.total == 3
    assert ocx._session.adapters["https://"].max_retries.allowed_methods is None


def test_model_classes(ocx, api_data):
    assert ocx.FunctionalProfiles


def test_sample_int_id(ocx, api_data):
    """
    Ensure that uuids do not get coerced into integers in the off chance that
    they're numeric-looking.

    Regression test for https://github.com/onecodex/onecodex/issues/200
    """
    sample = ocx.Samples.get("0111111111111110")
    assert sample.id == "0111111111111110"


def test_sample_get(ocx, api_data):
    sample = ocx.Samples.get("761bc54b97f64980")
    assert sample.size == 302369471
    assert sample.filename == "SRR2352223.fastq.gz"
    assert sample.__repr__() == '<Samples 761bc54b97f64980: "SRR2352223.fastq.gz">'
    assert isinstance(sample.created_at, datetime.datetime)

    classification = sample.primary_classification
    assert classification
    assert classification.complete

    tags = sample.tags
    assert len(tags) > 1
    assert "isolate" in [t.name for t in tags]


def test_sample_download(runner, ocx, api_data):
    with runner.isolated_filesystem():
        sample = ocx.Samples.get("761bc54b97f64980")
        sample.download()


def test_sample_download_awaiting_data(runner, ocx, api_data):
    with runner.isolated_filesystem():
        sample = ocx.Samples.get("761bc54b97f64980")
        sample.visibility = "awaiting data"
        with pytest.raises(OneCodexException):
            sample.download()


def test_document_download(runner, ocx, api_data):
    with runner.isolated_filesystem():
        doc = ocx.Documents.get("a4f6727a840a4df0")
        doc.download()


def test_download_without_filename(runner, ocx, api_data):
    with runner.isolated_filesystem():
        doc = ocx.Documents.get("a4f6727a840a4df0")

        with pytest.raises(OneCodexException, match="specify `path`, `file_obj`, or `_filename`"):
            doc._download("download_uri", _filename=None)


def test_download_path_exists(runner, ocx, api_data):
    with runner.isolated_filesystem():
        doc = ocx.Documents.get("a4f6727a840a4df0")
        doc.download()

        with pytest.raises(OneCodexException, match="already exists"):
            doc.download()


def test_download_use_potion_session(runner, ocx, api_data):
    with runner.isolated_filesystem():
        doc = ocx.Documents.get("a4f6727a840a4df0")
        doc._download("download_uri", doc.filename, use_potion_session=True)


def test_download_with_progressbar(runner, ocx, api_data):
    doc = ocx.Documents.get("a4f6727a840a4df0")

    with runner.isolated_filesystem():
        doc.download(progressbar=True)

    file_obj = io.BytesIO()
    doc.download(file_obj=file_obj, progressbar=True)


def test_download_file_obj(ocx, api_data):
    file_obj = io.BytesIO()
    doc = ocx.Documents.get("a4f6727a840a4df0")
    doc.download(file_obj=file_obj)

    file_obj.seek(0)
    data = file_obj.read()

    assert data == b'"1234567890"'


def test_resourcelist(ocx, api_data):
    sample = ocx.Samples.get("761bc54b97f64980")
    tags1 = onecodex.models.ResourceList(sample.tags._resource, onecodex.models.misc.Tags)

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
    assert "object of type" in str(e.value)

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
    assert "contain duplicate objects" in str(e.value)

    # SampleCollections can be added together
    new_samples = samples + other_samples
    assert len(samples) == len(other_samples) == 3
    assert len(new_samples) == 6

    # addition doesn't work with a SampleCollection and a lone Samples object
    single_sample = ocx.Samples.get("761bc54b97f64980")

    with pytest.raises(TypeError) as e:
        samples + single_sample
    assert "can only concatenate" in str(e.value)

    # you can make a new SampleCollection from a list of Samples
    onecodex.models.SampleCollection([s for s in all_samples[:7]])

    # and from a list of Classifications
    onecodex.models.SampleCollection([s.primary_classification for s in all_samples[:7]])

    # but not a combination of both
    with pytest.raises(OneCodexException) as e:
        onecodex.models.SampleCollection(
            [s.primary_classification for s in all_samples[:3]] + [s for s in all_samples[4:7]]
        )
    assert "but not both" in str(e.value)

    # and not using unwrapped potion resources
    with pytest.raises(OneCodexException) as e:
        onecodex.models.SampleCollection([s._resource for s in all_samples[:7]])
    assert "can only contain" in str(e.value)

    # test filtering of Samples in a collection
    new_collection = all_samples.filter(lambda s: s.filename.endswith("9.fastq.gz"))
    assert all([s.filename.endswith("9.fastq.gz") for s in new_collection])
    assert len(new_collection) == 7
    assert id(new_collection) != id(all_samples)


def test_dir_method(ocx, api_data):
    sample = ocx.Samples.get("761bc54b97f64980")

    instance_names = dir(sample)
    assert "where" not in instance_names  # we mask @classmethod's
    assert "created_at" in instance_names  # property on _resource
    assert "save" in instance_names  # function in py3, method in py2

    class_names = dir(ocx.Samples)
    assert "where" in class_names
    assert "created_at" not in class_names
    assert "save" in class_names  # instance methods are available off class


def test_get_failure_instructions(ocx):
    with pytest.raises(TypeError):
        ocx.Samples("direct_id")


def test_model_deletions(ocx, api_data):
    sample = ocx.Samples.get("761bc54b97f64980")
    sample.delete()

    analysis = sample.primary_classification
    with pytest.raises(MethodNotSupported):
        analysis.delete()


def test_model_updates(ocx, api_data):
    sample = ocx.Samples.get("761bc54b97f64980")
    sample.visibility = "public" if sample.visibility == "private" else "public"

    # Read-only field
    with pytest.raises(MethodNotSupported):
        sample.filename = "something_else"

    # No update resource
    analysis = sample.primary_classification
    with pytest.raises(MethodNotSupported):
        analysis.created_at = datetime.datetime.utcnow()


def test_sample_saving(ocx, api_data):
    sample = ocx.Samples.get("761bc54b97f64980")
    visibility = sample.visibility
    sample.visibility = "private" if visibility == "public" else "public"
    sample.save()
    assert sample.visibility is not visibility


def test_tag_saving_raises_exception(ocx, api_data):
    with mock.patch("onecodex.models.misc.Tags.where", return_value=[]):
        new_tag = ocx.Tags(name="new ta2")

    with pytest.raises(MethodNotSupported) as e:
        new_tag.save()

    assert "Tags cannot be saved directly" in str(e)


def test_metadata_saving(ocx, api_data):
    sample = ocx.Samples.get("761bc54b97f64980")
    metadata1 = sample.metadata
    metadata2 = ocx.Metadata.get("4fe05e748b5a4f0e")
    assert metadata1 == metadata2
    metadata1.description = "my new description -- testing!"
    metadata1.date_collected = datetime.datetime.now()
    metadata1.save()
    assert isinstance(metadata1.date_collected, datetime.datetime)
    assert metadata1.description == "my new description -- testing!"
    assert hasattr(metadata1, "sample")  # This will fail because we don't mock it in the return


def test_metadata_saving_should_ignore_updated_at_prop(ocx, api_data):
    m = ocx.Metadata.get("4fe05e748b5a4f0e")
    m.custom["foobar"] = "baz"
    m.save()
    patch_request = api_data.calls[-1].request
    assert patch_request.method == "PATCH"
    patch_data = json.loads(patch_request.body)
    assert "custom" in patch_data
    assert "updated_at" not in patch_data


def test_dir_patching(ocx, api_data):
    sample = ocx.Samples.get("761bc54b97f64980")
    props = {
        "id",
        "created_at",
        "filename",
        "visibility",
        "metadata",
        "owner",
        "primary_classification",
        "project",
        "size",
        "tags",
    }
    for prop in props:
        assert prop in dir(sample)
    assert len(sample.__dict__) == 1  # I'm not sure we *want* this...


def test_classification_methods(ocx, api_data):
    classification = ocx.Classifications.get("45a573fb7833449a")
    assert isinstance(classification, onecodex.models.analysis.Classifications)


# Sorting and where clauses
@pytest.mark.parametrize(
    "where_args,where_kwargs,queries",
    [
        ([], {"visibility": "public"}, ['where={"visibility": "public"}']),
        ([], {"visibility": "private"}, ['where={"visibility": "private"}']),
        ([], {"filename": "SRR1234.fastq.gz"}, ['where={"filename": "SRR1234.fastq.gz"}']),
        (
            [],
            {"filename": "SRR1234.fastq.gz", "sort": "visibility"},
            ['where={"filename": "SRR1234.fastq.gz"}'],
        ),
        (
            [],
            {"visibility": "private", "filename": "tmp.fa"},
            ['"filename": "tmp.fa"', '"visibility": "private"'],
        ),
        (
            [{"visibility": "private", "filename": "tmp.fa"}],
            {},
            ['"filename": "tmp.fa"', '"visibility": "private"'],
        ),
        (["761bc54b97f64980"], {}, ['"$uri": {"$in": ["/api/v1/samples/761bc54b97f64980"']),
        (
            ["/api/v1/samples/761bc54b97f64980"],
            {},
            ['"$uri": {"$in": ["/api/v1/samples/761bc54b97f64980"'],
        ),
    ],
)
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
        samples = ocx.Samples.search_public(filename="tmp.fa")
        assert len(samples) == 0
    samples = ocx.Samples.where(filename="tmp.fa", public=True)
    assert len(samples) == 0


def test_where_organization(ocx, api_data):
    samples = ocx.Samples.where(organization=True)
    assert len(samples) == 1


def test_public_project(ocx, api_data):
    with pytest.warns(DeprecationWarning):
        projs = ocx.Projects.search_public(name="One Codex Project")
        assert len(projs) == 0
    projs = ocx.Projects.where(name="One Codex Project", public=True)
    assert len(projs) == 0


def test_where_clauses_with_tags(ocx, api_data):
    tag = ocx.Tags.get("5c1e9e41043e4435")
    sample = ocx.Samples.get("7428cca4a3a04a8e")
    samples = ocx.Samples.where(tags=[tag])
    assert sample in samples

    query = '{"tags": {"$containsall": [{"$ref": "/api/v1/tags/5c1e9e41043e4435"}]}}'
    query_in_urls = []
    for c in responses.calls:
        url = unquote_plus(c.request.url)
        query_in_urls.append(query in url)

    assert any(query_in_urls)


def test_where_filter(ocx, api_data):
    samples = ocx.Samples.where(filter=lambda s: s.filename.endswith("9.fastq.gz"))
    assert all([s.filename.endswith("9.fastq.gz") for s in samples])
    assert len(samples) == 7

    # filter kwarg must be callable
    with pytest.raises(OneCodexException) as e:
        ocx.Samples.where(filter="not_callable")
    assert "Expected callable" in str(e.value)


def test_where_primary_classification(ocx, api_data):
    analysis = ocx.Analyses.get("935c2a3611944e39")
    sample = ocx.Samples.get("7428cca4a3a04a8e")
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
    assert a.sample.filename == "MSA-1000.16S.example.fastq.gz"
    assert a.job.name == "One Codex Database"
    assert a.sample.visibility == "public"


def test_jobs(ocx, api_data):
    jobs = ocx.Jobs.all()
    assert len(jobs) == 24

    jobs = ocx.Jobs.where(public=True)
    assert len(jobs) == 24
