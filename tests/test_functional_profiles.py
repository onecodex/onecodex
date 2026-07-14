import json

import pandas as pd

import pytest

from onecodex.models import SampleCollection
from onecodex.models import FunctionalProfiles


def test_query_for_functional_analysis(ocx, api_data):
    sample_id = "73b8349a30b04957"
    profile = ocx.FunctionalProfiles.where(sample=sample_id)
    assert isinstance(profile[0], FunctionalProfiles)

    profiles = ocx.FunctionalProfiles.all()
    assert len(profiles) == 3


def test_functional_profiles_table(ocx, api_data):
    func_profile = ocx.FunctionalProfiles.get("31ddae978aff475f")
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
    assert eggnog_df["taxon_name"].isna().all()

    all_df = func_profile.table(taxa_stratified=False)
    assert len(all_df) == 358
    assert all_df["taxon_name"].isna().all()
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


def test_functional_profiles_results(ocx, api_data):
    func_profile = ocx.FunctionalProfiles.get("eec4ac90d9104d1e")
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


def test_filtered_table_includes_taxon_columns_when_stratified(ocx, api_data):
    func_profile = ocx.FunctionalProfiles.get("31ddae978aff475f")

    stratified = func_profile.filtered_table(annotation="go", metric="rpk", taxa_stratified=True)
    assert not stratified.empty
    assert {"taxon_id", "taxon_name"} <= set(stratified.columns)  # is subset

    unstratified = func_profile.filtered_table(annotation="go", metric="rpk", taxa_stratified=False)
    assert not unstratified.empty
    assert "taxon_id" not in unstratified.columns
    assert "taxon_name" not in unstratified.columns


def test_filtered_table_empty_result_includes_taxon_columns(ocx, api_data, monkeypatch):
    func_profile = ocx.FunctionalProfiles.get("31ddae978aff475f")
    monkeypatch.setattr(
        func_profile,
        "_filtered_results",
        lambda **kwargs: {"table": [], "n_reads": 0, "n_mapped": 0},
    )

    stratified = func_profile.filtered_table(annotation="go", metric="rpk", taxa_stratified=True)
    assert stratified.empty
    # An empty table should still include "taxon_id" and "taxon_name" for downstream processing
    assert {"id", "name", "value", "taxon_id", "taxon_name"} <= set(stratified.columns)

    unstratified = func_profile.filtered_table(annotation="go", metric="rpk", taxa_stratified=False)
    assert unstratified.empty
    assert "taxon_id" not in unstratified.columns
    assert "taxon_name" not in unstratified.columns


def test_functional_profiles_fetch(ocx, api_data):
    sample_ids = [
        "543c9c046e3e4e09",
        "66c1531cb0b244f6",
        "37e5151e7bcb4f87",
        "7428cca4a3a04a8e",  # does not have functional results
    ]
    samples = [ocx.Samples.get(sample_id) for sample_id in sample_ids]
    sc = SampleCollection(samples, skip_missing=True)

    with pytest.warns(UserWarning, match="Functional profile not found.*7428cca4a3a04a8e"):
        # SampleCollection._functional_profiles_fetch() populates the .functional_profiles attribute cache
        functional_profiles = sc._functional_profiles

    assert len(functional_profiles) == 3
    for profile in functional_profiles:
        assert isinstance(profile, FunctionalProfiles)


def test_collate_functional_results(ocx, api_data):
    sample_ids = ["543c9c046e3e4e09", "66c1531cb0b244f6", "37e5151e7bcb4f87"]
    samples = [ocx.Samples.get(sample_id) for sample_id in sample_ids]
    sc = SampleCollection(samples)
    df, mapping = sc._functional_results(
        annotation="go", metric="rpk", taxa_stratified=True, fill_missing=False, filler=0
    )

    assert df.columns.name == "functional_profile_id"
    assert df.index.names == ["feature_id", "taxon_id", "taxon_name"]

    assert isinstance(df, pd.DataFrame)
    assert df.shape == (39, 3)
    assert len(mapping) == len(df.index.get_level_values("feature_id").unique())
    assert set(mapping.keys()) == set(df.index.get_level_values("feature_id"))

    df, mapping = sc._functional_results(
        annotation="eggnog", metric="cpm", taxa_stratified=False, fill_missing=True, filler=0
    )

    df, mapping = sc._functional_results(
        annotation="pathways", metric="coverage", taxa_stratified=True, fill_missing=False, filler=0
    )

    assert df.shape == (27, 3)
    assert len(mapping) == len(df.index.get_level_values("feature_id").unique())
    assert set(mapping.keys()) == set(df.index.get_level_values("feature_id"))

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

    result, _ = sc._functional_results(
        annotation="pfam", metric="cpm", taxa_stratified=False, fill_missing=False, filler=0
    )
    assert result.shape == (2, 3)
    result, _ = sc._functional_results(
        annotation="pfam", metric="rpk", taxa_stratified=False, fill_missing=False, filler=0
    )
    assert result.shape == (2, 3)
    result, _ = sc._functional_results(
        annotation="go", metric="rpk", taxa_stratified=True, fill_missing=False, filler=0
    )
    assert result.shape == (39, 3)


def test_to_df_for_functional_profiles(ocx, api_data):
    sample_ids = ["543c9c046e3e4e09", "66c1531cb0b244f6", "37e5151e7bcb4f87"]
    samples = [ocx.Samples.get(sample_id) for sample_id in sample_ids]
    sc = SampleCollection(samples)
    df = sc.to_df(analysis_type="functional")
    assert df.shape == (27, 3)
    df = sc.to_df(
        analysis_type="functional",
        annotation="eggnog",
        metric="cpm",
        taxa_stratified=False,
        fill_missing=True,
        filler=0,
    )
    assert df.shape == (7, 3)
    assert df.ocx_functional_group == "eggnog"
    assert df.ocx_metric == "cpm"
    assert df.ocx_metadata.shape == (3, 92)
    assert df.columns.name == "functional_profile_id"
    assert df.index.name == "feature_id"
    assert set(df.ocx_metadata["sample_id"]) == set(sample_ids)
    assert set(df.ocx_feature_name_map.keys()) == set(df.index)

    # Functional df doesn't have classification df attributes
    with pytest.raises(AttributeError):
        df.ocx_taxonomy
    with pytest.raises(AttributeError):
        df.ocx_rank

    with pytest.raises(ValueError):
        sc.to_df(analysis_type="foo")


def test_filter_functional_runs_to_newest_job(ocx, raw_api_data, custom_mock_requests):
    # Newer run
    raw_api_data["GET::api/v1/functional_profiles"] += [
        {
            "$uri": "/api/v1/functional_profiles/eec4ac90d9104d1f",
            "complete": True,
            "created_at": "2023-09-25T17:27:30.622286-07:00",
            "error_msg": "",
            "job": {"$ref": "/api/v1/jobs/59e7904ea8ed4244"},
            "sample": {"$ref": "/api/v1/samples/37e5151e7bcb4f87"},
            "cost": None,
            "dependencies": [],
            "draft": False,
            "success": True,
        }
    ]
    # Default job
    raw_api_data["GET::api/v1/jobs/59e7904ea8ed4202"] = {
        "$uri": "/api/v1/jobs/59e7904ea8ed4202",
        "analysis_type": "functional",
        "created_at": "2016-05-05T17:27:02.116480+00:00",
        "name": "Functional v1",
        "public": True,
    }
    # Newer job
    raw_api_data["GET::api/v1/jobs/59e7904ea8ed4244"] = {
        "$uri": "/api/v1/jobs/59e7904ea8ed4244",
        "analysis_type": "functional",
        "created_at": "2023-09-27T17:27:02.116480+00:00",
        "name": "Functional v2",
        "public": True,
    }

    with open("tests/data/api/v1/functional_profiles/bde18eb9407d4c2f/results/index.json") as fin:
        results = json.load(fin)
    raw_api_data[
        "GET::api/v1/functional_profiles/eec4ac90d9104d1f/filtered_results\\?functional_group=pathways&metric=coverage&taxa_stratified=True"
    ] = results

    with custom_mock_requests(raw_api_data):
        sample_ids = ["543c9c046e3e4e09", "66c1531cb0b244f6", "37e5151e7bcb4f87"]
        samples = [ocx.Samples.get(sample_id) for sample_id in sample_ids]
        sc = SampleCollection(samples)
        with pytest.warns(UserWarning, match="mixing functional profile versions"):
            df = sc.to_df(analysis_type="functional")

        # All samples are included (one column per functional profile)
        assert df.shape[1] == 3
        assert "eec4ac90d9104d1f" in df.columns
        # The newer version has correct PF00005 value for Campylobacter hominis
        key = ("PF00005", "76517", "g__Campylobacter.s__Campylobacter_hominis")
        assert df.loc[key, "eec4ac90d9104d1f"] == 256.524
