from hashlib import sha256
import pytest

pytest.importorskip("pandas")  # noqa
import warnings

from onecodex.exceptions import OneCodexException
from onecodex.lib.enums import Rank


def test_auto_rank(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    # auto rank should be species for shotgun data
    samples._collate_results(metric="abundance")
    assert samples._get_auto_rank("auto") == "species"

    samples._collate_results(metric="abundance_w_children")
    assert samples._get_auto_rank("auto") == "species"

    samples._collate_results(metric="readcount")
    assert samples._get_auto_rank("auto") == "species"

    samples._collate_results(metric="readcount_w_children")
    assert samples._get_auto_rank("auto") == "species"

    # inside the pandas extension, auto rank should choose that of the ResultsDataFrame
    results = samples.to_df(rank="phylum")
    assert results.ocx._get_auto_rank("auto") == "phylum"


def test_default_metric(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")
    assert samples._metric == "abundance_w_children"
    assert samples.metric == "Relative Abundance"


def test_all_nan_classification_id_property(ocx, api_data):
    import numpy as np

    samples = ocx.Samples.where(project="4b53797444f846c4")
    # should be the same when called on `SampleCollection` and on `OneCodexAccessor`
    assert samples._all_nan_classification_ids == []
    assert samples.to_df(fill_missing=False).ocx._all_nan_classification_ids == []

    samples._results.iloc[0] = np.nan

    assert samples.to_df(fill_missing=False).ocx._all_nan_classification_ids == [
        samples._results.iloc[0].name
    ]
    assert samples._all_nan_classification_ids == [samples._results.iloc[0].name]


def test_guess_normalization(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    norm_results = samples.to_df(normalize=True)
    assert norm_results.ocx._guess_normalized() is True

    with pytest.raises(OneCodexException) as e:
        unnorm_results = samples.to_df(normalize=False)
    assert "Data has already been normalized" in str(e.value)

    samples._collate_results(metric="readcount_w_children")
    assert samples._guess_normalized() is False

    norm_results = samples.to_df(normalize=True)
    assert norm_results.ocx._guess_normalized() is True

    unnorm_results = samples.to_df(normalize=False)
    assert unnorm_results.ocx._guess_normalized() is False


def test_normalization_with_zero_abundance_samples(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")
    samples._collate_results(metric="readcount_w_children")
    samples._results.iloc[:2, :] = 0
    assert not samples._guess_normalized()

    df = samples.to_df(normalize=True)

    assert list(df.sum(axis=1, skipna=False).round(6)) == [0.0, 0.0, 1.0]


def test_fill_missing_df(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")
    samples._collate_results(metric="abundance_w_children")

    assert samples.to_df(fill_missing=False).isnull().values.any() == True  # noqa
    assert samples.to_df(fill_missing=True).isnull().values.any() == False  # noqa


def test_metadata_fetch(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    samples._collate_results(metric="readcount_w_children")

    # tuple with invalid field
    with pytest.raises(OneCodexException) as e:
        samples._metadata_fetch([("does_not_exist", "vegetables")])
    assert "not found" in str(e.value)

    # tuple where one field is not categorical
    with pytest.raises(OneCodexException) as e:
        samples._metadata_fetch([("totalige", "vegetables")])
    assert "must be categorical" in str(e.value)

    # proper tuple of categorical metadata fields
    df, fields = samples._metadata_fetch([("eggs", "vegetables")])
    assert fields[("eggs", "vegetables")] == "eggs_vegetables"
    assert df["eggs_vegetables"].tolist() == ["True_True", "True_True", "True_True"]

    # Label reserved field name
    df, fields = samples._metadata_fetch(["Label"])
    assert fields["Label"] == "Label"
    assert df["Label"].tolist() == ["SRR4408293.fastq", "SRR4305031.fastq", "SRR4408292.fastq"]

    # exact match of single metadata field
    df, fields = samples._metadata_fetch(["totalige"])
    assert fields["totalige"] == "totalige"
    assert df["totalige"].tolist() == [62.9, 91.5, 112.0]

    # single metadata field doesn't match anything
    df, fields = samples._metadata_fetch(["does_not_exist"])
    assert fields["does_not_exist"] == "does_not_exist"
    assert df["does_not_exist"].tolist() == [None, None, None]

    # tax_id coerced to string from integer
    df, fields = samples._metadata_fetch([1279])
    assert fields[1279] == "Staphylococcus (1279)"
    assert df["Staphylococcus (1279)"].round(10).tolist() == [4.0172e-06, 4.89491e-05, 2.0881e-06]

    # tax_name, should take lowest matching tax_id, be case-insensitive, and report within-rank abundance
    df, fields = samples._metadata_fetch(["bacteroid"])
    assert fields["bacteroid"] == "Bacteroidaceae (815)"
    assert df["Bacteroidaceae (815)"].round(10).tolist() == [
        0.344903898,
        0.1656058794,
        0.7776093433,
    ]

    # label is a metadata field or callable
    df, fields = samples._metadata_fetch(["Label"], label="eggs")
    assert df["Label"].tolist() == ["True (1)", "True (2)", "True (3)"]

    df, fields = samples._metadata_fetch(["Label"], label=lambda x: str(x["eggs"]) + "_foo")
    assert df["Label"].tolist() == ["True_foo (1)", "True_foo (2)", "True_foo (3)"]

    # label must be a string or callable
    with pytest.raises(OneCodexException) as e:
        df, fields = samples._metadata_fetch(["Label"], label=False)
    assert "string or callable" in str(e.value)

    # label is a field not in the table
    with pytest.raises(OneCodexException) as e:
        df, fields = samples._metadata_fetch(["Label"], label="does_not_exist")
    assert "not found" in str(e.value)

    # label is a callable doesn't return a string
    with pytest.raises(OneCodexException) as e:
        df, fields = samples._metadata_fetch(["Label"], label=lambda x: False)
    assert "string from label" in str(e.value)


def test_results_filtering_rank(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    samples._collate_results(metric="readcount_w_children")

    tree = samples.tree_build()
    tree = samples.tree_prune_tax_ids(tree, ["1279"])

    # test filtering by rank
    genus_results = samples.to_df(rank="genus")
    family_results = samples.to_df(rank="family")

    assert "1279" in genus_results
    assert "1279" not in family_results
    assert "90964" in family_results
    assert "90964" not in genus_results

    # the resulting taxonomy dicts should have taxa at rank specified and their parents,
    # but no children, and no higher rank nodes not connected to taxa at the specified rank
    assert len(genus_results.ocx.taxonomy.loc[genus_results.ocx.taxonomy["rank"] == "species"]) == 0
    assert len(genus_results.ocx.taxonomy.loc[genus_results.ocx.taxonomy["rank"] == "genus"]) == 474
    assert (
        len(genus_results.ocx.taxonomy.loc[genus_results.ocx.taxonomy["rank"] == "family"]) == 163
    )
    assert (
        len(family_results.ocx.taxonomy.loc[family_results.ocx.taxonomy["rank"] == "species"]) == 0
    )
    assert len(family_results.ocx.taxonomy.loc[family_results.ocx.taxonomy["rank"] == "genus"]) == 0
    assert (
        len(family_results.ocx.taxonomy.loc[family_results.ocx.taxonomy["rank"] == "family"]) == 177
    )

    # catch invalid ranks
    with pytest.raises(OneCodexException, match="does_not_exist"):
        samples.to_df(rank="does_not_exist")

    # should warn if using rank=kingdom
    with warnings.catch_warnings(record=True) as w:
        samples.to_df(rank="kingdom")
        assert len(w) == 1
        assert "superkingdom" in str(w[-1].message)


def test_results_filtering_other(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    samples._collate_results(metric="readcount_w_children")

    # normalize the data
    norm_results = samples.to_df(normalize=True)
    assert norm_results.sum(axis=1).round(6).tolist() == [1.0, 1.0, 1.0]

    # trying to 'unnormalize' data fails
    with pytest.raises(OneCodexException) as e:
        norm_results.ocx.to_df(normalize=False)
    assert "already been normalized" in str(e.value)

    # remove columns where every value is zero
    results = samples.to_df(rank=None, normalize=False, remove_zeros=False)
    assert len(results.columns) == 3156
    results["1279"] = 0
    results["1280"] = 0
    results = results.ocx.to_df(rank=None, normalize=False, remove_zeros=True)
    assert len(results.columns) == 3154

    # return only taxa with at least 100 reads in one or more samples
    assert (
        (samples.to_df(rank=None, normalize=False, remove_zeros=False, threshold=100) >= 100).any(
            axis=0
        )
    ).all()

    # top N most abundant taxa
    assert (
        sha256(samples.to_df(top_n=10, rank="genus").round(6).to_json().encode()).hexdigest()
        == "437fcf282b885440571f552bf322056b6633154b247a9382ca074eee4b1ebf59"
    )

    # check entire contents of long and wide format tables
    wide_tbl = samples.to_df(table_format="wide", rank="genus", normalize=False)
    assert wide_tbl.sum().sum() == 38034443
    assert sum(wide_tbl.columns.astype(int).tolist()) == 109170075

    long_tbl = samples.to_df(table_format="long", rank="genus", normalize=False)
    assert long_tbl["Reads"].sum() == 38034443
    assert long_tbl["tax_id"].astype(int).sum() / 3 == 109170075

    # should fail if format is not wide or long
    with pytest.raises(OneCodexException) as e:
        samples.to_df(table_format="does_not_exist")
    assert "must be one of" in str(e.value)


@pytest.mark.parametrize("metric", ["readcount_w_children", "abundance_w_children"])
@pytest.mark.parametrize("rank", [Rank.Phylum, Rank.Genus])
def test_to_df_include_taxa_missing_rank(ocx, api_data, metric, rank):
    samples = ocx.Samples.where(project="4b53797444f846c4")
    samples._collate_results(metric=metric)

    df = samples.to_df(rank=rank, include_taxa_missing_rank=True)
    assert f"No {rank.value}" in df.columns


@pytest.mark.parametrize("metric", ["readcount", "abundance"])
def test_to_df_include_taxa_missing_rank_invalid_usage(ocx, api_data, metric):
    samples = ocx.Samples.where(project="4b53797444f846c4")
    samples._collate_results(metric=metric)

    with pytest.raises(OneCodexException, match="`rank`.*include_taxa_missing_rank"):
        samples.to_df(rank=None, include_taxa_missing_rank=True)

    with pytest.raises(OneCodexException, match="`include_taxa_missing_rank`.*metrics"):
        samples.to_df(include_taxa_missing_rank=True)
