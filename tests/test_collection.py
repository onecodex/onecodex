import json
from hashlib import sha256

import mock
import pytest

pytest.importorskip("pandas")  # noqa

from onecodex.exceptions import OneCodexException
from onecodex.lib.enums import Metric, Rank
from onecodex.models.collection import SampleCollection


def test_sample_collection_hashability(samples):
    # the entire sample collection is hashable
    assert len(set([samples])) == 1

    # individual samples cannot be hashed as they hold state (e.g., two samples with the same id
    # can have different metadata if one has been edited)
    with pytest.raises(TypeError):
        set(samples)


def test_sample_collection_pandas(samples):
    # manipulations of samples in the collection ought to update the stored dfs
    class_id = samples[2].primary_classification.id
    del samples[2]

    df = samples.to_df()
    assert df.ocx_metric is Metric.AbundanceWChildren

    assert len(samples) == 2
    assert len(df) == 2
    assert len(samples.metadata) == 2

    results, _ = samples._collate_results()

    assert class_id not in results.index
    assert class_id not in samples.metadata.index

    df = samples.to_df(metric=Metric.ReadcountWChildren)
    assert df.ocx_metric is Metric.ReadcountWChildren

    df = samples.to_df(metric=Metric.Readcount)
    assert df.ocx_metric is Metric.Readcount


def test_normalize_results(samples):
    """
    These numbers should match the numbers displayed on the
    classification results page (filtered=True): https://app.onecodex.com/analysis/6579e99943f84ad2
    """
    run = samples[0].primary_classification

    assert run.id == "6579e99943f84ad2"  # update golden output below if this ID changes

    assert run._classification_stats == {
        "n_host_reads": 109585,
        "n_mapped_microbial_reads": 24527027,
        "n_mapped_reads": 25444333,
        "n_nonspecific_reads": 807721,
        "n_reads_total": 32742710,
    }

    collection = SampleCollection([samples[0]])

    df_abundance = collection.to_df(rank=Rank.Species, metric=Metric.AbundanceWChildren)
    vals = df_abundance.values[0]
    vals.sort()
    assert vals[::-1][:6].round(4).tolist() == [0.2439, 0.0685, 0.0646, 0.0324, 0.0324, 0.0287]
    assert vals.sum().round(4) == 1.0

    df_readcount = collection.to_df(rank=Rank.Species, metric=Metric.Readcount)
    vals = df_readcount.values[0]
    vals.sort()
    assert vals[::-1][:6].round(4).tolist() == [4417519, 283510, 212457, 179055, 172214, 128158]

    df_readcount_w_children = collection.to_df(rank=Rank.Species, metric=Metric.ReadcountWChildren)
    vals = df_readcount_w_children.values[0]
    vals.sort()
    assert vals[::-1][:6].round(4).tolist() == [4417519, 754966, 283510, 226485, 214160, 212457]

    df_prop_classified_w_children = collection.to_df(
        rank=Rank.Species, metric=Metric.PropClassifiedWChildren
    )
    vals = df_prop_classified_w_children.values[0]
    vals.sort()
    assert vals[::-1][:6].round(4).tolist() == [0.1801, 0.0308, 0.0116, 0.0092, 0.0087, 0.0087]

    df_prop_classified_w_children = collection.to_df(
        rank=Rank.Genus, metric=Metric.PropClassifiedWChildren
    )
    vals = df_prop_classified_w_children.values[0]
    vals.sort()
    assert vals[::-1][:6].round(4).tolist() == [0.2424, 0.1801, 0.0404, 0.0396, 0.0313, 0.03]

    df_prop_classified_w_children = collection.to_df(
        rank=Rank.Phylum, metric=Metric.PropClassifiedWChildren
    )
    vals = df_prop_classified_w_children.values[0]
    vals.sort()
    assert vals[::-1][:6].round(4).tolist() == [0.5110, 0.2645, 0.1887, 0.012, 0.0008, 0.0000]

    df_prop_readcount = collection.to_df(rank=Rank.Species, metric=Metric.PropReadcount)
    vals = df_prop_readcount.values[0]
    vals.sort()
    assert vals[::-1][:6].round(4).tolist() == [0.1349, 0.0087, 0.0065, 0.0055, 0.0053, 0.0039]

    # not displayed on frontend, but consistent with data in Complete Results Table
    df_prop_readcount = collection.to_df(rank=Rank.Species, metric=Metric.PropReadcountWChildren)
    vals = df_prop_readcount.values[0]
    vals.sort()
    assert vals[::-1][:6].round(4).tolist() == [0.1349, 0.0231, 0.0087, 0.0069, 0.0065, 0.0065]


def test_automatic_rank(samples):
    assert samples.is_metagenomic
    assert samples.automatic_rank(metric=Metric.Auto) == "species"
    assert samples.automatic_rank(metric=Metric.AbundanceWChildren) == "species"
    assert samples.automatic_rank(metric=Metric.Abundance) == "species"
    assert samples.automatic_rank(metric=Metric.ReadcountWChildren) == "species"
    assert samples.automatic_rank(metric=Metric.Readcount) == "species"

    with mock.patch.object(SampleCollection, "is_metagenomic", property(lambda _: False)):
        assert samples.automatic_rank(metric=Metric.Auto) == "genus"
        assert samples.automatic_rank(metric=Metric.ReadcountWChildren) == "genus"
        assert samples.automatic_rank(metric=Metric.Readcount) == "genus"


@pytest.mark.filterwarnings("ignore:.*multiple analysis types.*")
def test_biom(ocx, api_data):
    c1 = ocx.Classifications.get("45a573fb7833449a")
    c2 = ocx.Classifications.get("593601a797914cbf")
    biom = SampleCollection([c1, c2]).to_otu(metric=Metric.ReadcountWChildren)
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

    # golden output
    assert biom["data"] == [[0, 0, 3], [0, 1, 80], [1, 0, 0], [1, 1, 0]]

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


def test_classifications(ocx, samples):
    # should work with a list of classifications as input, not just samples
    samples._oc_model = ocx.Classifications
    samples._res_list = samples._classifications
    assert len(samples._classifications) == len(samples)

    # should issue a warning if a classification did not succeed
    object.__setattr__(samples._res_list[0], "success", False)
    with pytest.warns(UserWarning, match="not successful"):
        samples._classifications


def test_classification_fetch_sample_missing_primary_classification(ocx, samples):
    sample = samples._res_list[0]
    object.__setattr__(sample, "primary_classification", None)

    msg = f"Classification not found.*{sample.id}"
    with pytest.warns(UserWarning, match=msg):
        samples._classifications

    samples._kwargs["skip_missing"] = False
    with pytest.raises(OneCodexException, match=msg):
        samples._classifications


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
        == "3ead672171efcb806323a55216683834aa89b5a657da31ab5bf01c6adcd882e6"
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
    results, taxonomy = samples._collate_results(metric=metric)

    # check contents of results df
    string_to_hash = ""
    for col in sorted(results.columns.tolist()):
        for row in sorted(results.index.tolist()):
            try:
                string_to_hash += results.fillna(0).loc[row, col].astype(str)
            except AttributeError:
                pass

    assert sha256(string_to_hash.encode()).hexdigest() == sha

    # check contents of taxonomy df
    string_to_hash = ""
    for col in sorted(taxonomy.columns.tolist()):
        for row in sorted(taxonomy.index.tolist()):
            try:
                string_to_hash += taxonomy.loc[row, col].astype(str)
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


def test_automatic_metric_no_abundance_estimates(samples_without_abundances):
    assert samples_without_abundances.automatic_metric == Metric.NormalizedReadcountWChildren


def test_automatic_metric_with_abundances(samples):
    assert samples.automatic_metric == Metric.AbundanceWChildren


def test_automatic_metric_majority_rules(samples, samples_without_abundances):
    assert (samples + samples_without_abundances[:1]).automatic_metric == Metric.AbundanceWChildren
    assert (
        samples_without_abundances + samples[:1]
    ).automatic_metric == Metric.NormalizedReadcountWChildren


def test_automatic_metric_majority_lacking_abundance_estimates(samples, samples_without_abundances):
    assert samples_without_abundances.automatic_metric == Metric.NormalizedReadcountWChildren
