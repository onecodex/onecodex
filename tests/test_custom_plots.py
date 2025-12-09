import io
import json
import uuid

import numpy as np
import pandas as pd
import pytest

from onecodex.exceptions import ValidationError
from onecodex.lib.enums import Metric
from onecodex.viz._custom_plots.collection import SampleCollection, Samples
from onecodex.viz._custom_plots.enums import ExportFormat, PlotRepr, PlotType
from onecodex.viz._custom_plots.metadata import _get_metadata_field_value
from onecodex.viz._custom_plots.models import PlotParams
from onecodex.viz._custom_plots.utils import get_plot_title


def generate_id() -> str:
    return uuid.uuid4().hex[0:16]


def load_classification_results_json(classification_uuid: str) -> dict | list:
    filepath = f"tests/data/api/v1/classifications/{classification_uuid}/results/index.json"
    with open(filepath, "r") as f:
        return json.load(f)


@pytest.fixture
def sample_collection() -> SampleCollection:
    sample1_uuid = generate_id()
    classification1_uuid = "0f4ee4ecb3a3412f"
    sample1 = Samples(
        {
            "uuid": sample1_uuid,
            "metadata": {
                "sample_id": sample1_uuid,
                "metadata_id": generate_id(),
                "classification_id": classification1_uuid,
                "filename": "sample1.fastq",
                "sample_name": "Sample 1",
                "sample_type": "Isolate",
                "library_type": "WGS",
                "platform": "Illumina",
                "cohort": "C1",
                "whitespace_field": " value 1",
            },
            "primary_classification": {
                "uuid": classification1_uuid,
                "job_name": "OCX DB",
                "api_results": load_classification_results_json(classification1_uuid),
            },
            "functional_profile": None,
        }
    )

    sample2_uuid = generate_id()
    classification2_uuid = "bef0bc57dd7f4c43"
    sample2 = Samples(
        {
            "uuid": sample2_uuid,
            "metadata": {
                "sample_id": sample2_uuid,
                "metadata_id": generate_id(),
                "classification_id": classification2_uuid,
                "filename": "sample2.fastq",
                "sample_name": "Sample 2",
                "sample_type": "Other",
                "library_type": "WGS",
                "platform": "Illumina",
                "cohort": "C2",
                "whitespace_field": "value 2  ",
            },
            "primary_classification": {
                "uuid": classification2_uuid,
                "job_name": "OCX DB",
                "api_results": load_classification_results_json(classification2_uuid),
            },
            "functional_profile": None,
        }
    )

    sample3_uuid = generate_id()
    classification3_uuid = "fadd35dd60074a5c"
    sample3 = Samples(
        {
            "uuid": sample3_uuid,
            "metadata": {
                "sample_id": sample3_uuid,
                "metadata_id": generate_id(),
                "classification_id": classification3_uuid,
                "filename": "sample3.fastq",
                "sample_name": "Sample 3",
                "sample_type": "Metagenomic",
                "library_type": "Other",
                "platform": "Other",
                "cohort": "C1",
                "whitespace_field": "  \tvalue 3  ",
            },
            "primary_classification": {
                "uuid": classification3_uuid,
                "job_name": "OCX DB",
                "api_results": load_classification_results_json(classification3_uuid),
            },
            "functional_profile": None,
        }
    )

    return SampleCollection([sample1, sample2, sample3])


@pytest.fixture
def default_plot_params_payload() -> dict:
    return {
        "tag": None,
        "project": None,
        "source_name": "Source",
        "plot_type": "taxa",
        "plot_repr": "bargraph",
        "metric": "abundance",
        "alpha_metric": "shannon",
        "beta_metric": "braycurtis",
        "export_format": None,
        "top_n": 10,
        "rank": "genus",
        "facet_by": None,
        "group_by": None,
        "secondary_group_by": None,
        "filter_by": None,
        "filter_value": [],
        "label_by": [],
        "sort_by": None,
        "functional_annotation": "pathways",
        "functional_metric": "rpk",
        "functional_pathways_metric": "rpk",
        "functional_top_n": 10,
        "functional_label": "name",
    }


@pytest.mark.parametrize(
    "params",
    [
        {"plot_type": "taxa", "plot_repr": "bargraph"},
        {
            "plot_type": "taxa",
            "plot_repr": "bargraph",
            "filter_by": "cohort",
            "filter_value": ["C1"],
        },
        {"plot_type": "taxa", "plot_repr": "heatmap"},
        {"plot_type": "alpha"},
        {"plot_type": "alpha", "group_by": "cohort", "secondary_group_by": "sample_type"},
        {"plot_type": "beta", "plot_repr": "pca"},
        {"plot_type": "beta", "plot_repr": "pcoa"},
        {"plot_type": "beta", "plot_repr": "distance"},
        {"metric": "auto"},
        {"export_format": None},
        {"export_format": "csv"},
        {"export_format": "xlsx"},
    ],
)
def test_plot(sample_collection, default_plot_params_payload, params):
    params = PlotParams.model_validate(default_plot_params_payload | params)
    result = sample_collection.plot(params)

    if params.metric == Metric.Auto:
        assert result.params != params
        assert result.params.metric == Metric.ReadcountWChildren
    else:
        assert result.params == params

    assert isinstance(result.chart, dict)
    assert isinstance(result.x_axis_label_links, dict)
    assert result.error is None
    assert result.warnings == []

    if params.export_format:
        assert isinstance(result.exported_chart_data, bytes)

        data = io.BytesIO(result.exported_chart_data)
        if params.export_format == ExportFormat.Csv:
            df = pd.read_csv(data)
        else:
            df = pd.read_excel(data)

        assert isinstance(df, pd.DataFrame)
    else:
        assert result.exported_chart_data is None


# Regression test for DEV-10319
#
# If a metadata name matches a taxonomic name and samples where the value for that metadata column
# all have NaN as the value, the metadata column would get matched to the taxonomy table causing a
# random taxon to get mapped to the plotting parameter, which then caused an error
def test_plot_does_not_match_on_taxonomy(sample_collection, default_plot_params_payload):
    sample_collection.taxonomy.loc["2", "name"] = "Group"
    sample_collection.taxonomy.loc["2", "rank"] = "An Invalid Rank"
    sample_collection.metadata["Group"] = "Test"
    sample_collection.metadata["Cohort"] = "Some Cohort"

    # These will get filtered out *after* the metadata key check
    sample_collection.metadata.loc[sample_collection.metadata.index[-3:], "Cohort"] = (
        "Some Other Cohort"
    )
    sample_collection.metadata.loc[sample_collection.metadata.index[-3:], "Group"] = np.nan

    params = PlotParams.model_validate(
        default_plot_params_payload
        | {
            "plot_type": "beta",
            "plot_repr": "pca",
            "facet_by": "Group",
            "filter_by": "Cohort",
            "filter_value": ["Some Other Cohort"],
            "rank": "species",
            "metric": "abundance_w_children",
        }
    )

    result = sample_collection.plot(params)
    assert result.error is None


@pytest.mark.parametrize(
    "params,columns",
    [
        (
            {"plot_type": "taxa", "plot_repr": "bargraph", "metric": "readcount_w_children"},
            {"Tax Name", "Sample 1", "Sample 2", "Sample 3"},
        ),
        (
            {"plot_type": "taxa", "plot_repr": "heatmap", "metric": "abundance"},
            {"Tax Name", "Sample 1", "Sample 2", "Sample 3"},
        ),
        (
            {
                "plot_type": "taxa",
                "plot_repr": "bargraph",
                "metric": "abundance",
                "group_by": "library_type",
            },
            {"Tax Name", "WGS", "Other"},
        ),
        (
            {
                "plot_type": "alpha",
                "metric": "readcount",
                "alpha_metric": "shannon",
            },
            {"Sample", "Classification ID", "Shannon"},
        ),
        (
            {
                "plot_type": "alpha",
                "metric": "abundance",
                "alpha_metric": "observed_taxa",
            },
            {"Sample", "Classification ID", "Observed Taxa"},
        ),
        (
            {
                "plot_type": "alpha",
                "group_by": "cohort",
                "metric": "readcount",
                "alpha_metric": "shannon",
            },
            {"Sample", "cohort", "Shannon"},
        ),
        (
            {
                "plot_type": "alpha",
                "facet_by": "cohort",
                "group_by": "platform",
                "metric": "readcount",
                "alpha_metric": "shannon",
            },
            {"Sample", "platform", "cohort", "Shannon"},
        ),
        (
            {
                "plot_type": "alpha",
                "group_by": "platform",
                "secondary_group_by": "cohort",
                "metric": "readcount",
                "alpha_metric": "shannon",
            },
            {"Sample", "platform", "cohort", "Shannon"},
        ),
        (
            {
                "plot_type": "alpha",
                "facet_by": "cohort",
                "group_by": "platform",
                "secondary_group_by": "cohort",
                "metric": "readcount",
                "alpha_metric": "shannon",
            },
            {"Sample", "platform", "cohort", "Shannon"},
        ),
        (
            {"plot_type": "beta", "plot_repr": "pca", "metric": "readcount_w_children"},
            {"Sample", "Classification ID", "PC1", "PC2"},
        ),
        (
            {"plot_type": "beta", "plot_repr": "pcoa", "metric": "abundance"},
            {"Sample", "Classification ID", "PC1", "PC2"},
        ),
        (
            {
                "plot_type": "beta",
                "plot_repr": "distance",
                "metric": "readcount_w_children",
            },
            {"Sample", "Sample 1", "Sample 2", "Sample 3"},
        ),
    ],
)
def test_export_chart_data(sample_collection, default_plot_params_payload, params, columns):
    params = PlotParams.model_validate(
        default_plot_params_payload | params | {"export_format": "csv"}
    )
    result = sample_collection.plot(params)
    df = pd.read_csv(io.BytesIO(result.exported_chart_data))

    assert set(df.columns) == columns

    if params.plot_repr == PlotRepr.Distance:
        assert set(df["Sample"]) == columns - {"Sample"}


def test_validate_plot_params_no_functional_profiles(
    sample_collection, default_plot_params_payload
):
    params = PlotParams.model_validate(default_plot_params_payload | {"plot_type": "functional"})

    with pytest.raises(ValidationError, match="Functional Analysis has not been run"):
        sample_collection._validate_plot_params(params)


@pytest.mark.parametrize(
    "params",
    [
        {"facet_by": "foo"},
        {"group_by": "foo"},
        {"secondary_group_by": "foo"},
        {"filter_by": "foo"},
        {"label_by": ["cohort", "foo"]},
        {"sort_by": "foo"},
    ],
)
def test_validate_plot_params_invalid_metadata_field(
    sample_collection, default_plot_params_payload, params
):
    params = PlotParams.model_validate(default_plot_params_payload | params)

    with pytest.raises(ValidationError, match="metadata field 'foo' does not exist"):
        sample_collection._validate_plot_params(params)


@pytest.mark.parametrize(
    "field,values,expected_num_samples",
    [
        ("cohort", ["C1"], 2),
        ("cohort", ["C1", "C2"], 3),
        ("sample_type", ["Metagenomic"], 1),
        ("sample_name", ["Sample 1", "Sample 3"], 2),
    ],
)
def test_filter_by_metadata(sample_collection, field, values, expected_num_samples):
    assert len(sample_collection) == 3

    filtered_collection = sample_collection._filter_by_metadata(field, values)

    assert len(filtered_collection) == expected_num_samples
    assert filtered_collection is not sample_collection


@pytest.mark.parametrize(
    "label_by,expected_labels",
    [
        (["cohort"], {"C1 (1)", "C1 (2)", "C2"}),
        (["cohort", "sample_type"], {"C1 / Isolate", "C1 / Metagenomic", "C2 / Other"}),
        (["whitespace_field"], {"value 1", "value 2", "value 3"}),
    ],
)
def test_x_axis_label_func(sample_collection, label_by, expected_labels):
    label_func = sample_collection._x_axis_label_func(PlotType.Taxa, label_by)

    labels = set()
    for record in sample_collection.metadata.to_dict("records"):
        labels.add(label_func(record))

    assert labels == expected_labels


def test_x_axis_label_classification_links(sample_collection):
    label_func = sample_collection._x_axis_label_func(PlotType.Taxa, ["cohort", "sample_type"])
    x_axis_label_links = sample_collection._x_axis_label_classification_links(label_func, None)

    for classification_id, record in sample_collection.metadata.to_dict("index").items():
        label = label_func(record)
        link = x_axis_label_links[label]
        assert link == f"/classification/{classification_id}"


def test_x_axis_label_classification_links_with_group_by(sample_collection):
    label_func = sample_collection._x_axis_label_func(PlotType.Taxa, ["cohort"])
    assert sample_collection._x_axis_label_classification_links(label_func, "group") == {}


def test_x_axis_sort_func(sample_collection):
    label_func = sample_collection._x_axis_label_func(PlotType.Taxa, ["cohort", "sample_type"])
    sort_func = sample_collection._x_axis_sort_func("sample_type", label_func)
    sorted_labels = sort_func(None)
    assert sorted_labels == ["C1 / Isolate", "C1 / Metagenomic", "C2 / Other"]


def test_x_axis_sort_func_no_sort_by(sample_collection):
    label_func = sample_collection._x_axis_label_func(PlotType.Taxa, ["sample_name"])
    sort_func = sample_collection._x_axis_sort_func(None, label_func)
    assert sort_func is None


@pytest.mark.parametrize(
    "metadata,field,expected",
    [
        ({}, "sample_name", "N/A"),
        ({}, "field1", "N/A"),
        ({"field1": np.nan}, "field1", "N/A"),
        ({"sample_name": "foo"}, "sample_name", "foo"),
        ({"field1": "foo"}, "field1", "foo"),
    ],
)
def test_get_metadata_field_value(metadata, field, expected):
    actual = _get_metadata_field_value(metadata, field)
    assert actual == expected


@pytest.mark.parametrize(
    "params,start",
    [
        ({"metric": "readcount"}, "Normalized readcount"),
        ({"metric": "readcount_w_children"}, "Normalized readcount with children"),
        ({"metric": "abundance"}, "Relative abundance"),
        ({"metric": "abundance_w_children"}, "Relative abundance"),
        ({"metric": "auto"}, "Taxa"),
        ({"facet_by": "Cohort"}, "Cohort"),
    ],
)
def test_get_plot_title(default_plot_params_payload, params, start):
    params = PlotParams.model_validate(
        default_plot_params_payload | params | {"source_name": "My Project"}
    )
    title = get_plot_title(params)
    assert title == f"{start} plot of My Project samples"
