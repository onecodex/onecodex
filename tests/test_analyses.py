from hashlib import sha256
import pytest

pytest.importorskip("pandas")  # noqa

from onecodex.exceptions import OneCodexException
from onecodex.lib.enums import Rank


def test_auto_rank(samples):
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


def test_default_metric(samples):
    assert samples._metric == "abundance_w_children"
    assert samples.metric == "Relative Abundance"


def test_classification_ids_without_abundances_property(samples):
    import numpy as np

    # should be the same when called on `SampleCollection` and on `OneCodexAccessor`
    assert samples._classification_ids_without_abundances == []
    assert samples.to_df(fill_missing=False).ocx._classification_ids_without_abundances == []

    samples._results.iloc[0] = np.nan

    assert samples.to_df(fill_missing=False).ocx._classification_ids_without_abundances == [
        samples._results.iloc[0].name
    ]
    assert samples._classification_ids_without_abundances == [samples._results.iloc[0].name]


def test_guess_normalization(samples):
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


def test_normalization_with_zero_abundance_samples(samples):
    samples._collate_results(metric="readcount_w_children")
    samples._results.iloc[:2, :] = 0
    assert not samples._guess_normalized()

    df = samples.to_df(normalize=True)

    assert list(df.sum(axis=1, skipna=False).round(6)) == [0.0, 0.0, 1.0]


def test_fill_missing_df(samples):
    samples._collate_results(metric="abundance_w_children")

    assert samples.to_df(fill_missing=False).isnull().values.any() == True  # noqa
    assert samples.to_df(fill_missing=True).isnull().values.any() == False  # noqa


def test_metadata_fetch(samples):
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
    results = samples._metadata_fetch([("eggs", "vegetables")])
    assert results.renamed_fields[("eggs", "vegetables")] == "eggs_vegetables"
    assert results.df["eggs_vegetables"].tolist() == ["True_True", "True_True", "True_True"]
    assert not results.taxonomy_fields

    # Label reserved field name
    results = samples._metadata_fetch(["Label"])
    assert results.renamed_fields["Label"] == "Label"
    assert results.df["Label"].tolist() == [
        "SRR4408293.fastq",
        "SRR4305031.fastq",
        "SRR4408292.fastq",
    ]
    assert not results.taxonomy_fields

    # exact match of single metadata field
    results = samples._metadata_fetch(["totalige"])
    assert results.renamed_fields["totalige"] == "totalige"
    assert results.df["totalige"].tolist() == [62.9, 91.5, 112.0]
    assert not results.taxonomy_fields

    # single metadata field doesn't match anything
    results = samples._metadata_fetch(["does_not_exist"])
    assert results.renamed_fields["does_not_exist"] == "does_not_exist"
    assert results.df["does_not_exist"].tolist() == [None, None, None]
    assert not results.taxonomy_fields

    # tax_id coerced to string from integer
    results = samples._metadata_fetch([1279])
    assert results.renamed_fields[1279] == "Staphylococcus (1279)"
    assert results.df["Staphylococcus (1279)"].round(10).tolist() == [
        4.0172e-06,
        4.89491e-05,
        2.0881e-06,
    ]
    assert results.taxonomy_fields == {1279}

    # tax_name, should take lowest matching tax_id, be case-insensitive, and report within-rank abundance
    results = samples._metadata_fetch(["bacteroid"])
    assert results.renamed_fields["bacteroid"] == "Bacteroidaceae (815)"
    assert results.df["Bacteroidaceae (815)"].round(10).tolist() == [
        0.344903898,
        0.1656058794,
        0.7776093433,
    ]
    assert results.taxonomy_fields == {"bacteroid"}

    # label is a metadata field or callable
    results = samples._metadata_fetch(["Label"], label="eggs")
    assert results.df["Label"].tolist() == ["True (1)", "True (2)", "True (3)"]
    assert not results.taxonomy_fields

    results = samples._metadata_fetch(["Label"], label=lambda x: str(x["eggs"]) + "_foo")
    assert results.df["Label"].tolist() == ["True_foo (1)", "True_foo (2)", "True_foo (3)"]
    assert not results.taxonomy_fields

    # label must be a string or callable
    with pytest.raises(OneCodexException) as e:
        samples._metadata_fetch(["Label"], label=False)
    assert "string or callable" in str(e.value)

    # label is a field not in the table
    with pytest.raises(OneCodexException) as e:
        samples._metadata_fetch(["Label"], label="does_not_exist")
    assert "not found" in str(e.value)

    # label is a callable doesn't return a string
    with pytest.raises(OneCodexException) as e:
        samples._metadata_fetch(["Label"], label=lambda x: False)
    assert "string from label" in str(e.value)


@pytest.mark.parametrize(
    "query",
    [
        1,  # exact match on tax ID, "no rank" rank
        "Bifidobacterium longum subsp. longum",  # exact match on taxon name, "subspecies" rank
        "cellular org",  # partial match on taxon name, "no rank" rank
    ],
)
def test_metadata_fetch_match_taxonomy_non_canonical_rank(samples, query):
    samples._collate_results(metric="readcount_w_children")
    results = samples._metadata_fetch([query], match_taxonomy=True)

    assert results.df[query].isnull().all()
    assert results.renamed_fields == {query: str(query)}
    assert results.taxonomy_fields == set()


def test_results_filtering_rank(samples):
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
    with pytest.warns(UserWarning, match="superkingdom"):
        samples.to_df(rank="kingdom")


def test_results_filtering_other(samples):
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
def test_to_df_include_taxa_missing_rank(samples, metric, rank):
    samples._collate_results(metric=metric)
    df = samples.to_df(rank=rank, include_taxa_missing_rank=True)
    assert f"No {rank.value}" in df.columns


@pytest.mark.parametrize("metric", ["readcount", "abundance"])
def test_to_df_include_taxa_missing_rank_invalid_usage(samples, metric):
    samples._collate_results(metric=metric)

    with pytest.raises(OneCodexException, match="`rank`.*include_taxa_missing_rank"):
        samples.to_df(rank=None, include_taxa_missing_rank=True)

    with pytest.raises(OneCodexException, match="`include_taxa_missing_rank`.*metrics"):
        samples.to_df(include_taxa_missing_rank=True)
