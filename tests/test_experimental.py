import pandas as pd

import pytest

from onecodex.models import SampleCollection
from onecodex.models.experimental import FunctionalProfiles


def test_query_for_functional_analysis(ocx_experimental, api_data):
    sample_id = "73b8349a30b04957"
    profile = ocx_experimental.FunctionalProfiles.where(sample=sample_id)
    assert isinstance(profile[0], FunctionalProfiles)


def test_functional_profiles_table(ocx_experimental, api_data):
    func_profile = ocx_experimental.FunctionalProfiles.get("31ddae978aff475f")
    df = func_profile.table()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 992
    assert set(df.columns) == {
        "group_name",
        "id",
        "metric",
        "name",
        "taxa_stratified",
        "taxon_id",
        "taxon_name",
        "value",
    }
    assert len(df["taxon_name"].unique()) == 47

    eggnog_df = func_profile.table(annotation="eggnog", taxa_stratified=False)
    assert set(eggnog_df["group_name"]) == {"eggnog"}
    assert list(eggnog_df["taxon_name"].unique()) == [None]

    all_df = func_profile.table(taxa_stratified=False)
    assert len(all_df) == 358
    assert list(all_df["taxon_name"].unique()) == [None]
    assert list(all_df["group_name"].unique()) == [
        "gene_family",
        "metacyc",
        "pfam",
        "pathways",
        "go",
        "eggnog",
        "reaction",
        "ec",
        "ko",
    ]


def test_functional_profiles_results(ocx_experimental, api_data):
    func_profile = ocx_experimental.FunctionalProfiles.get("eec4ac90d9104d1e")
    json_results = func_profile.results()
    assert set(json_results.keys()) == {"n_mapped", "n_reads", "table"}
    assert isinstance(json_results["table"], list)
    assert set(json_results["table"][0].keys()) == {
        "group_name",
        "id",
        "metric",
        "name",
        "taxa_stratified",
        "taxon_id",
        "taxon_name",
        "value",
    }


def test_functional_profiles_fetch(ocx_experimental, api_data):
    sample_ids = ["543c9c046e3e4e09", "66c1531cb0b244f6", "37e5151e7bcb4f87"]
    samples = [ocx_experimental.Samples.get(sample_id) for sample_id in sample_ids]
    sc = SampleCollection(samples)
    # SampleCollection._functional_profiles_fetch() populates the .functional_profiles attribute cache
    functional_profiles = sc._functional_profiles
    for profile in functional_profiles:
        assert isinstance(profile, FunctionalProfiles)


def test_collate_functional_results(ocx_experimental, api_data):
    sample_ids = ["543c9c046e3e4e09", "66c1531cb0b244f6", "37e5151e7bcb4f87"]
    samples = [ocx_experimental.Samples.get(sample_id) for sample_id in sample_ids]
    sc = SampleCollection(samples)
    df, mapping = sc._functional_results(
        annotation="go", metric="rpk", taxa_stratified=True, fill_missing=False, filler=0
    )
    assert isinstance(df, pd.DataFrame)
    assert df.shape == (3, 40)  # 39 features + functional profile id
    assert len(mapping) == 39
    # The first column is the functional profile id
    assert sorted(list(mapping.keys())) == sorted(list(df.columns)[1:])
    assert df.compare(
        sc._cached[
            "functional_results_annotation=go_metric=rpk_taxa_stratified=True_fill_missing=False_filler=0"
        ]
    ).empty
    df, mapping = sc._functional_results(
        annotation="eggnog", metric="cpm", taxa_stratified=False, fill_missing=True, filler=0
    )
    # Old cache is still kept
    assert (
        "functional_results_annotation=go_metric=rpk_taxa_stratified=True_fill_missing=False_filler=0"
        in sc._cached
    )
    assert df.compare(
        sc._cached[
            "functional_results_annotation=eggnog_metric=cpm_taxa_stratified=False_fill_missing=True_filler=0"
        ]
    ).empty
    df, mapping = sc._functional_results(
        annotation="pathways", metric="coverage", taxa_stratified=True, fill_missing=False, filler=0
    )
    assert df.shape == (3, 28)
    assert len(mapping) == 27
    # The first column is the functional profile id
    assert sorted(list(mapping.keys())) == sorted(list(df.columns)[1:])
    with pytest.raises(ValueError):
        sc._functional_results(
            annotation="all", metric="rpk", taxa_stratified=True, fill_missing=False, filler=0
        )
    with pytest.raises(ValueError):
        sc._functional_results(
            annotation="pathways", metric="rpk", taxa_stratified=True, fill_missing=False, filler=0
        )
    with pytest.raises(ValueError):
        sc._functional_results(
            annotation="go", metric="coverage", taxa_stratified=False, fill_missing=True, filler=0
        )
    sc._functional_results(
        annotation="pfam", metric="cpm", taxa_stratified=False, fill_missing=False, filler=0
    )
    assert sc._cached[
        "functional_results_annotation=pfam_metric=cpm_taxa_stratified=False_fill_missing=False_filler=0"
    ].shape == (3, 3)
    sc._functional_results(
        annotation="pfam", metric="rpk", taxa_stratified=False, fill_missing=False, filler=0
    )
    assert sc._cached[
        "functional_results_annotation=pfam_metric=rpk_taxa_stratified=False_fill_missing=False_filler=0"
    ].shape == (3, 3)
    sc._functional_results(
        annotation="go", metric="rpk", taxa_stratified=True, fill_missing=False, filler=0
    )
    assert sc._cached[
        "functional_results_annotation=go_metric=rpk_taxa_stratified=True_fill_missing=False_filler=0"
    ].shape == (3, 40)


def test_to_df_for_functional_profiles(ocx_experimental, api_data):
    sample_ids = ["543c9c046e3e4e09", "66c1531cb0b244f6", "37e5151e7bcb4f87"]
    samples = [ocx_experimental.Samples.get(sample_id) for sample_id in sample_ids]
    sc = SampleCollection(samples)
    df = sc.to_df(analysis_type="functional")
    assert df.shape == (3, 28)
    df = sc.to_df(
        analysis_type="functional",
        annotation="eggnog",
        metric="cpm",
        taxa_stratified=False,
        fill_missing=True,
        filler=0,
    )
    assert df.shape == (3, 8)
    assert df.ocx_functional_group == "eggnog"
    assert df.ocx_metric == "cpm"
    assert df.ocx_metadata.shape == (3, 92)
    assert df.index.name == "sample_id"
    assert list(df.index.values) == sample_ids
    assert list(df.ocx_metadata["sample_id"]) == sample_ids
    # The first column is the functional profile id
    assert sorted(list(df.ocx_feature_name_map.keys())) == sorted(list(df.columns)[1:])

    # Functional df doesn't have classification df attributes
    with pytest.raises(AttributeError):
        df.ocx_taxonomy
    with pytest.raises(AttributeError):
        df.ocx_rank

    with pytest.raises(ValueError):
        sc.to_df(analysis_type="foo")
