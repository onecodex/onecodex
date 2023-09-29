import json

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
    assert df.shape == (3, 39)
    assert len(mapping) == 39
    assert sorted(list(mapping.keys())) == sorted(list(df.columns))
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
    assert df.shape == (3, 27)
    assert len(mapping) == 27
    assert sorted(list(mapping.keys())) == sorted(list(df.columns))
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
    ].shape == (3, 2)
    sc._functional_results(
        annotation="pfam", metric="rpk", taxa_stratified=False, fill_missing=False, filler=0
    )
    assert sc._cached[
        "functional_results_annotation=pfam_metric=rpk_taxa_stratified=False_fill_missing=False_filler=0"
    ].shape == (3, 2)
    sc._functional_results(
        annotation="go", metric="rpk", taxa_stratified=True, fill_missing=False, filler=0
    )
    assert sc._cached[
        "functional_results_annotation=go_metric=rpk_taxa_stratified=True_fill_missing=False_filler=0"
    ].shape == (3, 39)


def test_to_df_for_functional_profiles(ocx_experimental, api_data):
    sample_ids = ["543c9c046e3e4e09", "66c1531cb0b244f6", "37e5151e7bcb4f87"]
    samples = [ocx_experimental.Samples.get(sample_id) for sample_id in sample_ids]
    sc = SampleCollection(samples)
    df = sc.to_df(analysis_type="functional")
    assert df.shape == (3, 27)
    df = sc.to_df(
        analysis_type="functional",
        annotation="eggnog",
        metric="cpm",
        taxa_stratified=False,
        fill_missing=True,
        filler=0,
    )
    assert df.shape == (3, 7)
    assert df.ocx_functional_group == "eggnog"
    assert df.ocx_metric == "cpm"
    assert df.ocx_metadata.shape == (3, 92)
    assert df.index.name == "sample_id"
    assert set(df.index.values) == set(sample_ids)
    assert set(df.ocx_metadata["sample_id"]) == set(sample_ids)
    assert set(df.ocx_feature_name_map.keys()) == set(df.columns)

    # Functional df doesn't have classification df attributes
    with pytest.raises(AttributeError):
        df.ocx_taxonomy
    with pytest.raises(AttributeError):
        df.ocx_rank

    with pytest.raises(ValueError):
        sc.to_df(analysis_type="foo")


def test_filter_functional_runs_to_newest_job(ocx_experimental, raw_api_data, custom_mock_requests):
    for json_profile in raw_api_data["GET::api/v1_experimental/functional_profiles"]:
        json_profile["job"]["$ref"] = json_profile["job"]["$ref"].replace("v1/", "v1_experimental/")
    # Newer run
    raw_api_data["GET::api/v1_experimental/functional_profiles"] += [
        {
            "$uri": "/api/v1_experimental/functional_profiles/eec4ac90d9104d1f",
            "complete": True,
            "created_at": "2023-09-25T17:27:30.622286-07:00",
            "error_msg": "",
            "job": {"$ref": "/api/v1_experimental/jobs/59e7904ea8ed4244"},
            "sample": {"$ref": "/api/v1/samples/37e5151e7bcb4f87"},
            "success": True,
        }
    ]
    # Default job
    raw_api_data["GET::api/v1_experimental/jobs/59e7904ea8ed4202"] = {
        "$uri": "/api/v1_experimental/jobs/59e7904ea8ed4202",
        "analysis_type": "functional",
        "created_at": "2016-05-05T17:27:02.116480+00:00",
        "name": "Functional v1",
        "public": True,
    }
    # Newer job
    raw_api_data["GET::api/v1_experimental/jobs/59e7904ea8ed4244"] = {
        "$uri": "/api/v1_experimental/jobs/59e7904ea8ed4244",
        "analysis_type": "functional",
        "created_at": "2023-09-27T17:27:02.116480+00:00",
        "name": "Functional v2",
        "public": True,
    }

    with open(
        "tests/data/api/v1_experimental/functional_profiles/bde18eb9407d4c2f/results/index.json"
    ) as fin:
        results = json.load(fin)
    raw_api_data[
        "GET::api/v1_experimental/functional_profiles/eec4ac90d9104d1f/filtered_results\\?functional_group=%22pathways%22&metric=%22coverage%22&taxa_stratified=true"
    ] = results

    with custom_mock_requests(raw_api_data):
        sample_ids = ["543c9c046e3e4e09", "66c1531cb0b244f6", "37e5151e7bcb4f87"]
        samples = [ocx_experimental.Samples.get(sample_id) for sample_id in sample_ids]
        sc = SampleCollection(samples)
        df = sc.to_df(analysis_type="functional")
        # All samples are included, one with newer version has proper values
        assert df.shape == (3, 112)
        assert df.loc["37e5151e7bcb4f87", "PF00005"] == 256.524  # older version has 4919.47
