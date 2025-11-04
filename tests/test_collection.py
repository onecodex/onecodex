from hashlib import sha256
import json

import pytest

pytest.importorskip("pandas")  # noqa

from onecodex.models.collection import SampleCollection
from onecodex.exceptions import OneCodexException


def test_sample_collection_pandas(samples):
    # manipulations of samples in the collection ought to update the stored dfs
    class_id = samples[2].primary_classification.id
    del samples[2]

    assert len(samples) == 2
    assert len(samples.to_df()) == 2
    assert len(samples.metadata) == 2
    assert class_id not in samples._results.index
    assert class_id not in samples.metadata.index


def test_collection_constructor(samples):
    col = SampleCollection(samples)
    assert col._kwargs["metric"] == "auto"

    with pytest.deprecated_call():
        col = SampleCollection(samples, field="readcount_w_children")
    assert isinstance(col, SampleCollection)
    assert col._kwargs["metric"] == "readcount_w_children"

    col = SampleCollection(samples, metric="abundance_w_children")
    assert isinstance(col, SampleCollection)
    assert col._kwargs["metric"] == "abundance_w_children"

    with pytest.raises(OneCodexException):
        SampleCollection(samples, metric="abundance_w_children", field="readcount_w_children")


@pytest.mark.filterwarnings("ignore:.*multiple analysis types.*")
def test_biom(ocx, api_data):
    c1 = ocx.Classifications.get("45a573fb7833449a")
    c2 = ocx.Classifications.get("593601a797914cbf")
    biom = SampleCollection([c1, c2], ocx.Classifications).to_otu()
    assert set(biom.keys()) == {
        "columns",
        "data",
        "date",
        "format",
        "format_url",
        "generated_by",
        "id",
        "matrix_element_type",
        "matrix_type",
        "rows",
        "shape",
        "type",
    }
    assert biom["format_url"] == "http://biom-format.org"
    assert set(biom["columns"][0].keys()) == {"id", "sample_id", "sample_filename", "metadata"}
    assert set(biom["rows"][0].keys()) == {"id", "metadata"}
    assert "taxonomy" in biom["rows"][0]["metadata"]

    # IDs
    assert len(biom["columns"]) == 2

    assert biom["columns"][0]["id"] == c1.id
    assert biom["columns"][0]["sample_id"] == c1.sample.id

    # Reults
    assert len(biom["rows"]) == 2

    assert biom["rows"][0]["metadata"]["taxonomy"] == [
        "",
        "",
        "",
        "",
        "",
        "",
        "Staphylococcus",
        "Staphylococcus sp. HGB0015",
    ]
    assert biom["rows"][1]["metadata"]["taxonomy"] == [
        "",
        "",
        "",
        "",
        "",
        "",
        "Staphylococcus",
        "",
    ]

    # make sure data is the correct shape
    assert [len(x) for x in biom["data"]] == [3, 3, 3, 3]

    # Format is row_id, sample id (column), count (sparse)
    assert biom["data"][0] == [0, 0, 3]
    assert biom["data"][1] == [0, 1, 80]
    assert biom["data"][2] == [1, 0, 0]
    assert biom["data"][3] == [1, 1, 0]
    assert biom["rows"][0]["id"] == "1078083"
    assert biom["rows"][1]["id"] == "1279"

    # Check that we're not including unnecessary references
    assert "$uri" not in biom["columns"][0]["metadata"]
    assert "sample" not in biom["columns"][0]["metadata"]

    # Test serialization
    assert json.loads(json.dumps(biom)) == biom  # tests json serialization


def test_classification_fetch(ocx, samples):
    # should work with a list of classifications as input, not just samples
    samples._oc_model = ocx.Classifications
    samples._res_list = samples._classifications
    samples._classification_fetch()

    # should issue a warning if a classification did not succeed
    object.__setattr__(samples._res_list[0], "success", False)
    samples._cached = {}
    with pytest.warns(UserWarning, match="not successful"):
        samples._classification_fetch()


def test_classification_fetch_sample_missing_primary_classification(ocx, samples):
    sample = samples._res_list[0]
    object.__setattr__(sample, "primary_classification", None)

    msg = f"Classification not found.*{sample.id}"
    with pytest.warns(UserWarning, match=msg):
        samples._classification_fetch()

    samples._kwargs["skip_missing"] = False
    with pytest.raises(OneCodexException, match=msg):
        samples._classification_fetch()


def test_collate_metadata(samples):
    # check contents of metadata df--at least that which can easily be coerced to strings
    metadata = samples.metadata
    string_to_hash = ""
    for col in sorted(metadata.columns.tolist()):
        for row in sorted(metadata.index.tolist()):
            val = metadata.loc[row, col]
            if isinstance(val, str):
                string_to_hash += val
                continue

            try:
                string_to_hash += val.astype(str)
            except AttributeError:
                pass

    assert (
        sha256(string_to_hash.encode()).hexdigest()
        == "d886c5850325327050ef0d73f942c92f8830524363b7a5cbf1eaed805a9f4582"
    )


@pytest.mark.parametrize(
    "metric,sha",
    [
        (
            "readcount_w_children",
            "dbe3adf601ca9584a49b1b5fcb1873dec5ea33986afa3f614e96609f9320c8ba",
        ),
        (
            "abundance_w_children",
            "6e08480e867e21a15d7e36a3bae4f772d20b89c1937e6587865b29b66374c483",
        ),
    ],
)
def test_collate_results(samples, metric, sha):
    samples._collate_results(metric=metric)

    # check contents of results df
    string_to_hash = ""
    for col in sorted(samples._results.columns.tolist()):
        for row in sorted(samples._results.index.tolist()):
            try:
                string_to_hash += samples._results.fillna(0).loc[row, col].astype(str)
            except AttributeError:
                pass

    assert sha256(string_to_hash.encode()).hexdigest() == sha

    # check contents of taxonomy df
    string_to_hash = ""
    for col in sorted(samples.taxonomy.columns.tolist()):
        for row in sorted(samples.taxonomy.index.tolist()):
            try:
                string_to_hash += samples.taxonomy.loc[row, col].astype(str)
            except AttributeError:
                pass

    assert (
        sha256(string_to_hash.encode()).hexdigest()
        == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    )

    # invalid field name
    with pytest.raises(OneCodexException) as e:
        samples._collate_results(metric="does_not_exist")
    assert "not valid" in str(e.value)


def test_collate_results_with_auto_metric_no_abundance_estimates(samples):
    assert samples._metric == "abundance_w_children"
    classifications = samples._classifications

    for classification in classifications:
        results = classification.results()
        for data in results["table"]:
            data["abundance_w_children"] = None
            data["abundance"] = None

    samples._collate_results(metric="auto")
    assert samples._metric == "readcount_w_children"


def test_collate_results_with_auto_metric_one_lacking_abundance_estimates(samples):
    assert samples._metric == "abundance_w_children"
    classifications = samples._classifications

    c1_results = classifications[0].results()
    for data in c1_results["table"]:
        data["abundance_w_children"] = None
        data["abundance"] = None

    samples._collate_results(metric="auto")
    # only 1 out of 3 does not have abundance estimates so we can
    # still use `abundance_w_children`
    assert samples._metric == "abundance_w_children"


def test_collate_results_with_auto_metric_majority_lacking_abundance_estimates(samples):
    assert samples._metric == "abundance_w_children"
    classifications = samples._classifications

    for classification in classifications[:2]:
        results = classification.results()
        for data in results["table"]:
            data["abundance_w_children"] = None
            data["abundance"] = None

    samples._collate_results(metric="auto")
    assert samples._metric == "readcount_w_children"
