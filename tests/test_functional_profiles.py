import json
import os
from collections import Counter

import pandas as pd
import pytest
import requests

from onecodex.exceptions import OneCodexException
from onecodex.lib.enums import FunctionalAnnotations, FunctionalAnnotationsMetric
from onecodex.models import FunctionalProfiles, SampleCollection


def _normalize_functional_rows(rows, comparison_keys, annotation=None):
    normalized_rows = []

    for row in rows:
        normalized_row = row.copy()
        row_annotation = normalized_row.get("group_name", annotation)

        # Original API results combined MetaCyc IDs and names, whereas condensed
        # results separate them where possible.
        if row_annotation == "metacyc" and normalized_row["name"] is not None:
            normalized_row["id"] = f"{normalized_row['id']}: {normalized_row['name']}"
            normalized_row["name"] = None

        normalized_rows.append(tuple(normalized_row[key] for key in comparison_keys))

    return normalized_rows


def test_query_for_functional_analysis(ocx, api_data):
    sample_id = "73b8349a30b04957"
    profile = ocx.FunctionalProfiles.where(sample=sample_id)
    assert isinstance(profile[0], FunctionalProfiles)

    profiles = ocx.FunctionalProfiles.all()
    assert len(profiles) == 3


def test_functional_profiles_table(ocx, api_data):
    func_profile = ocx.FunctionalProfiles.get("31ddae978aff475f")
    df = func_profile.table()
    # no metric is passed to table then it should be all metrics
    assert set(df["metric"]) == {"cpm", "rpk", "coverage", "abundance"}
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 536
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
    # no metric is passed to table then it should be all available metrics for that annotation
    assert set(eggnog_df["metric"]) == {"cpm", "rpk"}
    assert set(eggnog_df["group_name"]) == {"eggnog"}
    assert eggnog_df["taxon_name"].isna().all()

    all_df = func_profile.table(taxa_stratified=False)
    assert len(all_df) == 128
    assert all_df["taxon_name"].isna().all()
    assert list(all_df["group_name"].unique()) == [
        "eggnog",
        "go",
        "ko",
        "ec",
        "pfam",
        "reaction",
        "metacyc",
        "pathways",
    ]


@pytest.mark.parametrize(
    ("annotation", "metric"),
    [
        ("eggnog", "rpk"),
        ("metacyc", "cpm"),
        ("pathways", "coverage"),
    ],
)
@pytest.mark.parametrize("taxa_stratified", [False, True])
def test_functional_profiles_table_filters(ocx, api_data, annotation, metric, taxa_stratified):
    func_profile = ocx.FunctionalProfiles.get("31ddae978aff475f")

    df = func_profile.table(
        annotation=annotation,
        metric=metric,
        taxa_stratified=taxa_stratified,
    )

    assert not df.empty
    assert set(df["group_name"]) == {annotation}
    assert set(df["metric"]) == {metric}
    assert set(df["taxa_stratified"]) == {taxa_stratified}


def test_functional_profiles_table_filters_metric_without_annotation(ocx, api_data):
    func_profile = ocx.FunctionalProfiles.get("31ddae978aff475f")

    df = func_profile.table(metric="cpm", taxa_stratified=False)

    assert not df.empty
    assert set(df["metric"]) == {"cpm"}
    assert set(df["group_name"]) == {
        "eggnog",
        "go",
        "ko",
        "ec",
        "pfam",
        "reaction",
        "metacyc",
    }


@pytest.mark.parametrize("taxa_stratified", [False, True])
def test_functional_profiles_table_complete_abundance(ocx, api_data, taxa_stratified):
    func_profile = ocx.FunctionalProfiles.get("a888fdc70221befa")
    condensed = func_profile._condensed_results()

    complete_pathway_ids = {
        pathway[0] for pathway in condensed["results"]["pathways"] if pathway[3] == 1.0
    }

    df = func_profile.table(
        annotation="pathways",
        metric="complete_abundance",
        taxa_stratified=taxa_stratified,
    )

    assert not df.empty
    assert set(df["group_name"]) == {"pathways"}
    assert set(df["metric"]) == {"abundance"}
    assert set(df["id"]) <= complete_pathway_ids
    assert set(df["taxa_stratified"]) == {taxa_stratified}


def test_functional_profiles_table_rejects_invalid_metric(ocx, api_data):
    func_profile = ocx.FunctionalProfiles.get("31ddae978aff475f")

    with pytest.raises(
        OneCodexException,
        match="metric cpm cannot be retrieved for functional group pathways",
    ):
        func_profile.table(annotation="pathways", metric="cpm")


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
        FunctionalProfiles,
        "table",
        lambda self, **kwargs: pd.DataFrame(
            columns=["id", "name", "value", "taxon_id", "taxon_name"]
        ),
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
            "results_uri": os.path.abspath(
                "tests/data/api/v1/functional_profiles/bde18eb9407d4c2f/results/results_api_data.v1.json.gz"
            ),
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
            df = sc.to_df(analysis_type="functional", annotation="pfam", metric="cpm")

        # All samples are included (one row per functional profile)
        assert df.shape[0] == 3
        assert "eec4ac90d9104d1f" in df.index
        # The newer version has correct PF00005 value for Campylobacter hominis (76517)
        key = ("PF00005", "76517")
        assert df.loc["eec4ac90d9104d1f", key] == 256.524


def test_rehydrate_condensed_functional_results(ocx, api_data, original_functional_api_results):
    profile = ocx.FunctionalProfiles.get("a888fdc70221befa")

    assert profile.results_uri is not None

    condensed = profile._condensed_results()
    assert condensed is not None

    actual = profile.results()
    expected = original_functional_api_results

    assert actual["n_reads"] == expected["n_reads"]
    assert actual["n_mapped"] == expected["n_mapped"]

    recoverable_keys = (
        "group_name",
        "id",
        "name",
        "metric",
        "value",
        "taxa_stratified",
        "taxon_id",
    )

    actual_rows = _normalize_functional_rows(actual["table"], recoverable_keys)
    expected_rows = _normalize_functional_rows(expected["table"], recoverable_keys)

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
    ocx, api_data, annotation, metric, taxa_stratified, original_functional_results_filtered
):
    profile = ocx.FunctionalProfiles.get("a888fdc70221befa")
    condensed = profile._condensed_results()

    assert condensed is not None

    original_groups = tuple(condensed["results"])

    request_count = len(api_data.calls)
    actual = profile.filtered_table(
        annotation=annotation,
        metric=metric,
        taxa_stratified=taxa_stratified,
    )
    # actual_table = actual.to_dict(orient="records")
    # the previous /filtered_results endpoint removed unmapped/ungrouped/unintegrated
    actual_table = [
        row
        for row in actual.astype(object).where(pd.notna(actual), None).to_dict(orient="records")
        if row["id"] not in {"UNMAPPED", "UNGROUPED", "UNINTEGRATED"}
    ]

    # condensed results should be loaded locally without calling the API endpoint
    assert len(api_data.calls) == request_count

    # Fetch the same selection from the mocked filtered-results API endpoint
    expected = original_functional_results_filtered(
        annotation=annotation,
        metric=metric,
        taxa_stratified=taxa_stratified,
    )

    expected_keys = {"id", "name", "value"}
    comparison_keys = ("id", "name", "value")

    if taxa_stratified:
        expected_keys |= {"taxon_id", "taxon_name"}
        comparison_keys += ("taxon_id",)

    assert all(set(row) == expected_keys for row in actual_table)
    assert all(set(row) == expected_keys for row in expected["table"])

    actual_rows = _normalize_functional_rows(actual_table, comparison_keys, annotation=annotation)
    expected_rows = _normalize_functional_rows(
        expected["table"], comparison_keys, annotation=annotation
    )

    assert Counter(actual_rows) == Counter(expected_rows)

    if taxa_stratified:
        taxa_map = {node["id"]: node.get("name") for node in condensed["taxonomy"]["nodes"]}
        taxa_map["0"] = "unclassified"

        for row in actual_table:
            assert row["taxon_name"] == taxa_map.get(row["taxon_id"])

    # Filtering must not alter the cached condensed payload.
    assert tuple(condensed["results"]) == original_groups


@pytest.mark.parametrize(
    ("annotation", "metric"),
    [
        (annotation, metric)
        for annotation in FunctionalAnnotations
        for metric in FunctionalAnnotationsMetric.metrics_for_annotation(annotation)
    ],
)
@pytest.mark.parametrize("taxa_stratified", [False, True])
def test_to_functional_df_with_condensed_results(
    ocx, api_data, monkeypatch, annotation, metric, taxa_stratified
):
    profile = ocx.FunctionalProfiles.get("a888fdc70221befa")
    sample = ocx.Samples.get("37e5151e7bcb4f87")
    collection = SampleCollection([sample])

    expected = profile.filtered_table(
        annotation=annotation,
        metric=metric,
        taxa_stratified=taxa_stratified,
    )
    # these are now included in the table results
    expected = expected[~expected["id"].isin({"UNMAPPED", "UNGROUPED", "UNINTEGRATED"})]

    # use the standalone condensed profile instead of the profile returned by the
    # existing sample fixture.
    collection.__dict__["_functional_profiles"] = [profile]

    df = collection.to_functional_df(
        annotation=annotation,
        metric=metric,
        taxa_stratified=taxa_stratified,
        fill_missing=False,
    )

    # expected_values = {(row["id"], row["taxon_id"]): row["value"] for row in expected["table"]}
    if taxa_stratified:
        keys = zip(expected["id"], expected["taxon_id"])
        assert df.columns.names == ["feature_id", "taxon_id"]
    else:
        keys = expected["id"]
        assert df.columns.name == "feature_id"

    expected_values = dict(zip(keys, expected["value"]))

    assert df.shape == (1, len(expected_values))
    assert list(df.index) == [profile.id]
    assert df.loc[profile.id].to_dict() == expected_values
    assert df.ocx_functional_group == annotation
    assert df.ocx_metric == metric


@pytest.mark.parametrize("failure_mode", ["expired_url", "unsupported_version"])
def test_condensed_results_dont_fall_back_to_api(
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

    assert profile._condensed_results() is None

    request_count = len(api_data.calls)

    result = profile.table(
        annotation="go",
        metric="rpk",
        taxa_stratified=False,
    )

    assert result.empty

    assert len(api_data.calls) == request_count
