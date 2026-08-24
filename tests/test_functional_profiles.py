import json
from collections import Counter

import pandas as pd
import pytest
import requests

from onecodex.models import FunctionalProfiles, SampleCollection


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

    assert df.index.name == "functional_profile_id"
    assert df.columns.names == ["feature_id", "taxon_id"]

    assert isinstance(df, pd.DataFrame)
    assert df.shape == (3, 39)
    assert len(mapping) == len(df.columns.get_level_values("feature_id").unique())
    assert set(mapping.keys()) == set(df.columns.get_level_values("feature_id"))

    df, mapping = sc._functional_results(
        annotation="eggnog", metric="cpm", taxa_stratified=False, fill_missing=True, filler=0
    )

    df, mapping = sc._functional_results(
        annotation="pathways", metric="coverage", taxa_stratified=True, fill_missing=False, filler=0
    )

    assert df.shape == (3, 27)
    assert len(mapping) == len(df.columns.get_level_values("feature_id").unique())
    assert set(mapping.keys()) == set(df.columns.get_level_values("feature_id"))

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
    assert result.shape == (3, 2)
    result, _ = sc._functional_results(
        annotation="pfam", metric="rpk", taxa_stratified=False, fill_missing=False, filler=0
    )
    assert result.shape == (3, 2)
    result, _ = sc._functional_results(
        annotation="go", metric="rpk", taxa_stratified=True, fill_missing=False, filler=0
    )
    assert result.shape == (3, 39)


def test_to_df_for_functional_profiles(ocx, api_data):
    sample_ids = ["543c9c046e3e4e09", "66c1531cb0b244f6", "37e5151e7bcb4f87"]
    samples = [ocx.Samples.get(sample_id) for sample_id in sample_ids]
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
    assert df.index.name == "functional_profile_id"
    assert df.columns.name == "feature_id"
    assert set(df.ocx_metadata["sample_id"]) == set(sample_ids)
    assert set(df.ocx_feature_name_map.keys()) == set(df.columns)

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
        "draft": False,
    }
    # Newer job
    raw_api_data["GET::api/v1/jobs/59e7904ea8ed4244"] = {
        "$uri": "/api/v1/jobs/59e7904ea8ed4244",
        "analysis_type": "functional",
        "created_at": "2023-09-27T17:27:02.116480+00:00",
        "name": "Functional v2",
        "public": True,
        "draft": False,
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

        # All samples are included (one row per functional profile)
        assert df.shape[0] == 3
        assert "eec4ac90d9104d1f" in df.index
        # The newer version has correct PF00005 value for Campylobacter hominis (76517)
        key = ("PF00005", "76517")
        assert df.loc["eec4ac90d9104d1f", key] == 256.524


def test_rehydrate_condensed_functional_results(ocx, api_data):
    profile = ocx.FunctionalProfiles.get("a888fdc70221befa")

    assert profile.results_uri is not None

    condensed = profile._condensed_results()
    assert condensed is not None

    actual = profile.results()

    # the results directly from the API endpoint
    response = profile._client.get(f"{profile._api._base_url}{profile.field_uri}/results")
    response.raise_for_status()
    expected = response.json()

    assert actual["n_reads"] == expected["n_reads"]
    assert actual["n_mapped"] == expected["n_mapped"]

    def normalize(rows):
        recoverable_keys = (
            "group_name",
            "id",
            "name",
            "metric",
            "value",
            "taxa_stratified",
            "taxon_id",
        )

        return [tuple(row[key] for key in recoverable_keys) for row in rows]

    actual_rows = normalize(actual["table"])
    expected_rows = normalize(expected["table"])

    # all rows in results['table'] must match exactly
    assert Counter(actual_rows) == Counter(expected_rows)

    # the order of all functional groups annotations should be the same except for metacyc
    assert [row for row in actual_rows if row[0] != "metacyc"] == [
        row for row in expected_rows if row[0] != "metacyc"
    ]

    taxa_map = {node["id"]: node.get("name") for node in condensed["taxonomy"]["nodes"]}
    taxa_map["0"] = "unclassified"

    # condensed results reconstructs taxa names from its own taxonomy rather than using
    # the clade strings which are no longer supported
    for row in actual["table"]:
        expected_taxon_name = taxa_map.get(row["taxon_id"]) if row["taxa_stratified"] else None
        assert row["taxon_name"] == expected_taxon_name


@pytest.mark.parametrize(
    ("annotation", "metric"),
    [
        ("pathways", "abundance"),
        ("pathways", "complete_abundance"),
        ("pathways", "coverage"),
        ("metacyc", "cpm"),
        ("metacyc", "rpk"),
        ("eggnog", "cpm"),
        ("eggnog", "rpk"),
        ("go", "cpm"),
        ("go", "rpk"),
        ("ko", "cpm"),
        ("ko", "rpk"),
        ("ec", "cpm"),
        ("ec", "rpk"),
        ("pfam", "cpm"),
        ("pfam", "rpk"),
        ("reaction", "cpm"),
        ("reaction", "rpk"),
    ],
)
@pytest.mark.parametrize(
    "taxa_stratified",
    [
        False,
        True,
    ],
)
def test_rehydrate_condensed_filtered_functional_results(
    ocx,
    api_data,
    annotation,
    metric,
    taxa_stratified,
):
    profile = ocx.FunctionalProfiles.get("a888fdc70221befa")
    condensed = profile._condensed_results()

    assert condensed is not None

    original_groups = tuple(condensed["results"])

    request_count = len(api_data.calls)
    actual = profile._filtered_results(
        annotation=annotation,
        metric=metric,
        taxa_stratified=taxa_stratified,
    )

    # condensed results should be loaded locally without calling the API endpoint
    assert len(api_data.calls) == request_count

    # Fetch the same selection from the mocked filtered-results API endpoint
    response = profile._client.get(
        f"{profile._api._base_url}{profile.field_uri}/filtered_results",
        params={
            "functional_group": annotation,
            "metric": metric,
            "taxa_stratified": taxa_stratified,
        },
    )
    response.raise_for_status()
    expected = response.json()

    expected_keys = {"id", "name", "value"}
    comparison_keys = ("id", "name", "value")

    if taxa_stratified:
        expected_keys |= {"taxon_id", "taxon_name"}
        comparison_keys += ("taxon_id",)

    assert all(set(row) == expected_keys for row in actual["table"])
    assert all(set(row) == expected_keys for row in expected["table"])

    def normalize(rows):
        return [tuple(row[key] for key in comparison_keys) for row in rows]

    assert Counter(normalize(actual["table"])) == Counter(normalize(expected["table"]))
    assert actual["n_reads"] == expected["n_reads"]
    assert actual["n_mapped"] == expected["n_mapped"]

    if taxa_stratified:
        taxa_map = {node["id"]: node.get("name") for node in condensed["taxonomy"]["nodes"]}
        taxa_map["0"] = "unclassified"

        for row in actual["table"]:
            assert row["taxon_name"] == taxa_map.get(row["taxon_id"])

    # Filtering must not alter the cached condensed payload.
    assert tuple(condensed["results"]) == original_groups


def test_to_functional_df_with_condensed_results(ocx, api_data):
    profile = ocx.FunctionalProfiles.get("a888fdc70221befa")
    sample = ocx.Samples.get("37e5151e7bcb4f87")

    collection = SampleCollection([sample])

    # use the standalone condensed profile instead of the profile returned by the
    # existing sample fixture.
    collection.__dict__["_functional_profiles"] = [profile]

    request_count = len(api_data.calls)

    df = collection.to_functional_df(
        annotation="pathways",
        metric="coverage",
        taxa_stratified=True,
        fill_missing=False,
    )

    # Metadata may be fetched, but functional results must not come from
    # the filtered-results endpoint.
    new_requests = api_data.calls[request_count:]
    assert all("/filtered_results" not in call.request.url for call in new_requests)

    response = profile._client.get(
        f"{profile._api._base_url}{profile.field_uri}/filtered_results",
        params={
            "functional_group": "pathways",
            "metric": "coverage",
            "taxa_stratified": True,
        },
    )
    response.raise_for_status()

    expected_values = {
        (row["id"], row["taxon_id"]): row["value"] for row in response.json()["table"]
    }

    assert df.shape == (1, 556)
    assert list(df.index) == [profile.id]
    assert df.columns.names == ["feature_id", "taxon_id"]
    assert df.loc[profile.id].to_dict() == expected_values
    assert df.ocx_functional_group == "pathways"
    assert df.ocx_metric == "coverage"


@pytest.mark.parametrize("failure_mode", ["expired_url", "unsupported_version"])
def test_condensed_results_fall_back_to_api(
    ocx,
    api_data,
    monkeypatch,
    failure_mode,
):
    profile = ocx.FunctionalProfiles.get("a888fdc70221befa")

    if failure_mode == "expired_url":

        def load_results(_):
            raise requests.HTTPError("expired results URL")

    else:

        def load_results(_):
            return {"version": 2}

    monkeypatch.setattr("onecodex.models.analysis._load_results_uri", load_results)

    # make sure nothing from a previous test is cached
    FunctionalProfiles._condensed_results.cache_clear()

    try:
        assert profile._condensed_results() is None

        request_count = len(api_data.calls)

        result = profile._filtered_results(
            annotation="go",
            metric="rpk",
            taxa_stratified=False,
        )
    finally:
        # just prevent any further caching
        FunctionalProfiles._condensed_results.cache_clear()

    assert len(api_data.calls) == request_count + 1
    # we should be falling back to the api for results
    assert "/filtered_results" in api_data.calls[-1].request.url
    assert len(result["table"]) == 2952
    assert result["n_reads"] == 5334942
    assert result["n_mapped"] == 4128346
