import gzip
import io
import json
import os
import uuid

import mock
import numpy as np
import pandas as pd
import pytest

from onecodex.exceptions import ValidationError
from onecodex.lib.enums import (
    AdjustmentMethod,
    AlphaDiversityMetric,
    AlphaDiversityStatsTest,
    BetaDiversityMetric,
    BetaDiversityStatsTest,
    FunctionalAnnotations,
    FunctionalAnnotationsMetric,
    FunctionalLabel,
    Metric,
    PosthocStatsTest,
    Rank,
)
from onecodex.stats import (
    AlphaDiversityStatsResults,
    AncombcResults,
    BetaDiversityStatsResults,
    PosthocResults,
)
from onecodex.viz._custom_plots.collection import SampleCollection, Samples
from onecodex.viz._custom_plots.enums import ExportFormat, PlotRepr, PlotType, StatsType
from onecodex.viz._custom_plots.metadata import _get_metadata_field_value
from onecodex.viz._custom_plots.models import PlotParams, PlotResults, StatsParams, StatsResults
from onecodex.viz._custom_plots.utils import get_plot_title


def generate_id() -> str:
    return uuid.uuid4().hex[0:16]


def load_classification_results_json(classification_uuid: str) -> dict | list:
    base = f"tests/data/api/v1/classifications/{classification_uuid}/raw_results/index.json"
    if os.path.exists(base + ".gz"):
        with gzip.open(base + ".gz", "rt") as f:
            return json.load(f)
    with open(base, "r") as f:
        return json.load(f)


JOB_UUID = "a1b2c3d4e5f6a7b8"


def make_sample(
    *,
    classification_uuid: str,
    job_uuid: str = JOB_UUID,
    job_name: str = "OCX DB",
    functional_profile: dict | None = None,
    **extra_metadata,
) -> Samples:
    sample_uuid = generate_id()
    classification_id = generate_id()
    return Samples(
        {
            "uuid": sample_uuid,
            "metadata": {
                "sample_id": sample_uuid,
                "metadata_id": generate_id(),
                "classification_id": classification_id,
                **extra_metadata,
            },
            "primary_classification": {
                "uuid": classification_id,
                "job_uuid": job_uuid,
                "job_name": job_name,
                "api_results": load_classification_results_json(classification_uuid),
            },
            "functional_profile": functional_profile,
        }
    )


@pytest.fixture
def sample_collection() -> SampleCollection:
    sample1 = make_sample(
        classification_uuid="0f4ee4ecb3a3412f",
        filename="sample1.fastq",
        sample_name="Sample 1",
        sample_type="Isolate",
        library_type="WGS",
        platform="Illumina",
        cohort="C1",
        Bacteroides=None,
        whitespace_field=" value 1",
    )

    sample2 = make_sample(
        classification_uuid="bef0bc57dd7f4c43",
        filename="sample2.fastq",
        sample_name="Sample 2",
        sample_type="Other",
        library_type="WGS",
        platform="Illumina",
        cohort="C2",
        Bacteroides=None,
        whitespace_field="value 2\t",
    )

    sample3 = make_sample(
        classification_uuid="fadd35dd60074a5c",
        filename="sample3.fastq",
        sample_name="Sample 3",
        sample_type="Metagenomic",
        library_type="Other",
        platform="Other",
        cohort="C1",
        Bacteroides=None,
        whitespace_field="  \tvalue 3  ",
    )

    return SampleCollection([sample1, sample2, sample3])


@pytest.fixture
def sample_collection_mixed_abundances(sample_collection) -> SampleCollection:
    """A collection mixing samples with and without abundance estimates."""
    no_abundance_sample = make_sample(
        classification_uuid="464a7ebcnocaffe1",
        filename="no_abundances.fastq",
        sample_name="No Abundances",
        cohort="C1",
    )
    return SampleCollection([*sample_collection.samples, no_abundance_sample])


@pytest.fixture
def default_plot_params_payload() -> dict:
    return {
        "tag": None,
        "project": None,
        "source_name": "Source",
        "plot_type": PlotType.Taxa,
        "plot_repr": PlotRepr.Bargraph,
        "metric": Metric.Auto,
        "alpha_metric": AlphaDiversityMetric.Shannon,
        "beta_metric": BetaDiversityMetric.BrayCurtis,
        "export_format": None,
        "top_n": 10,
        "rank": Rank.Auto,
        "facet_by": None,
        "group_by": None,
        "secondary_group_by": None,
        "filter_by": None,
        "filter_value": [],
        "label_by": [],
        "sort_by": None,
        "functional_annotation": FunctionalAnnotations.Pathways,
        "functional_metric": FunctionalAnnotationsMetric.Rpk,
        "functional_pathways_metric": FunctionalAnnotationsMetric.Rpk,
        "functional_top_n": 10,
        "functional_label": FunctionalLabel.Name,
    }


@pytest.mark.parametrize(
    "params",
    [
        {"plot_type": PlotType.Taxa, "plot_repr": PlotRepr.Bargraph},
        {
            "plot_type": PlotType.Taxa,
            "plot_repr": PlotRepr.Bargraph,
            "filter_by": "cohort",
            "filter_value": ["C1"],
        },
        {
            "plot_type": PlotType.Taxa,
            "plot_repr": PlotRepr.Heatmap,
            "metric": Metric.NormalizedReadcountWChildren,
        },
        {"plot_type": PlotType.Alpha, "metric": Metric.NormalizedReadcountWChildren},
        {
            "plot_type": PlotType.Alpha,
            "group_by": "cohort",
            "secondary_group_by": "sample_type",
            "metric": Metric.NormalizedReadcountWChildren,
        },
        {
            "plot_type": PlotType.Beta,
            "plot_repr": PlotRepr.Pca,
            "metric": Metric.NormalizedReadcountWChildren,
        },
        {"plot_type": PlotType.Beta, "plot_repr": PlotRepr.Pcoa},
        {"plot_type": PlotType.Beta, "plot_repr": PlotRepr.Distance},
        {"metric": Metric.Auto},
        {"export_format": None},
        {"export_format": ExportFormat.Csv},
        {"export_format": ExportFormat.Xlsx},
    ],
)
def test_plot(sample_collection, default_plot_params_payload, params):
    params = PlotParams.model_validate(default_plot_params_payload | params)

    result = sample_collection.plot(params)

    assert result.error is None

    if params.metric == Metric.Auto:
        assert result.params.metric == sample_collection.automatic_metric
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


@pytest.mark.parametrize(
    "metric",
    [
        Metric.RawReadcount,
        Metric.RawReadcountWChildren,
        Metric.NormalizedRawReadcount,
        Metric.NormalizedRawReadcountWChildren,
        Metric.Abundance,
        Metric.AbundanceWChildren,
    ],
)
def test_plot_no_warning_metrics(
    sample_collection_mixed_abundances, default_plot_params_payload, metric
):
    """Raw readcount and abundance metrics must not emit a PlottingWarning for mixed-abundance
    collections."""
    params = PlotParams.model_validate(default_plot_params_payload | {"metric": metric})
    result = sample_collection_mixed_abundances.plot(params)

    assert result.error is None
    assert result.warnings == [], f"Unexpected warnings for metric {metric}: {result.warnings}"


@pytest.mark.parametrize(
    "metric",
    [
        Metric.Readcount,
        Metric.ReadcountWChildren,
        Metric.NormalizedReadcount,
        Metric.NormalizedReadcountWChildren,
    ],
)
def test_plot_warns_for_filtered_metrics_with_mixed_abundances(
    sample_collection_mixed_abundances, default_plot_params_payload, metric
):
    """Filtered readcount metrics should emit a PlottingWarning for mixed-abundance collections."""
    params = PlotParams.model_validate(default_plot_params_payload | {"metric": metric})
    result = sample_collection_mixed_abundances.plot(params)

    assert result.error is None
    assert any(
        "no abundances calculated" in w for w in result.warnings
    ), f"Expected warning for metric {metric}, got: {result.warnings}"


# Regression test for DEV-10319
#
# If a metadata name matches a taxonomic name and samples where the value for that metadata column
# all have NaN as the value, the metadata column would get matched to the taxonomy table causing a
# random taxon to get mapped to the plotting parameter, which then caused an error
def test_plot_does_not_match_on_taxonomy(sample_collection, default_plot_params_payload):
    params = PlotParams.model_validate(
        default_plot_params_payload
        | {
            "plot_type": PlotType.Beta,
            "plot_repr": PlotRepr.Pca,
            "facet_by": "Bacteroides",  # all are None
            "filter_by": "cohort",
            "filter_value": ["C1", "C2"],
            "rank": Rank.Genus,
            "metric": Metric.ReadcountWChildren,
        }
    )

    df = sample_collection.to_df(rank=Rank.Genus, metric=Metric.ReadcountWChildren)

    # make sure this taxon is actually there
    assert (sample_collection.taxonomy["name"] == "Bacteroides").sum() == 1
    assert "816" in df.columns

    result = sample_collection.plot(params)
    assert result.error is None


@pytest.mark.parametrize(
    "params,columns",
    [
        *[
            (
                {"plot_type": PlotType.Taxa, "plot_repr": PlotRepr.Bargraph, "metric": metric},
                {"Tax Name", "Sample 1", "Sample 2", "Sample 3"},
            )
            for metric in [
                Metric.Readcount,
                Metric.ReadcountWChildren,
                Metric.NormalizedReadcount,
                Metric.NormalizedReadcountWChildren,
                Metric.Abundance,
                Metric.AbundanceWChildren,
            ]
        ],
        *[
            (
                {
                    "plot_type": PlotType.Taxa,
                    "plot_repr": PlotRepr.Bargraph,
                    "metric": metric,
                    "group_by": "library_type",
                },
                {"Tax Name", "WGS", "Other"},
            )
            for metric in [
                Metric.Readcount,
                Metric.ReadcountWChildren,
                Metric.NormalizedReadcount,
                Metric.NormalizedReadcountWChildren,
                Metric.Abundance,
            ]
        ],
        (
            {"plot_type": PlotType.Taxa, "plot_repr": PlotRepr.Heatmap, "metric": Metric.Abundance},
            {"Tax Name", "Sample 1", "Sample 2", "Sample 3"},
        ),
        (
            {
                "plot_type": PlotType.Alpha,
                "metric": Metric.Readcount,
                "alpha_metric": AlphaDiversityMetric.Shannon,
            },
            {"Sample", "Classification ID", "Shannon"},
        ),
        (
            {
                "plot_type": PlotType.Alpha,
                "metric": Metric.Abundance,
                "alpha_metric": AlphaDiversityMetric.ObservedTaxa,
            },
            {"Sample", "Classification ID", "Observed Taxa"},
        ),
        (
            {
                "plot_type": PlotType.Alpha,
                "group_by": "cohort",
                "metric": Metric.Readcount,
                "alpha_metric": AlphaDiversityMetric.Shannon,
            },
            {"Sample", "cohort", "Shannon"},
        ),
        (
            {
                "plot_type": PlotType.Alpha,
                "facet_by": "cohort",
                "group_by": "platform",
                "metric": Metric.Readcount,
                "alpha_metric": AlphaDiversityMetric.Shannon,
            },
            {"Sample", "platform", "cohort", "Shannon"},
        ),
        (
            {
                "plot_type": PlotType.Alpha,
                "group_by": "platform",
                "secondary_group_by": "cohort",
                "metric": Metric.Readcount,
                "alpha_metric": AlphaDiversityMetric.Shannon,
            },
            {"Sample", "platform", "cohort", "Shannon"},
        ),
        (
            {
                "plot_type": PlotType.Alpha,
                "facet_by": "cohort",
                "group_by": "platform",
                "secondary_group_by": "cohort",
                "metric": Metric.Readcount,
                "alpha_metric": AlphaDiversityMetric.Shannon,
            },
            {"Sample", "platform", "cohort", "Shannon"},
        ),
        (
            {
                "plot_type": PlotType.Beta,
                "plot_repr": PlotRepr.Pca,
                "metric": Metric.ReadcountWChildren,
            },
            {"Sample", "Classification ID", "PC1", "PC2"},
        ),
        (
            {"plot_type": PlotType.Beta, "plot_repr": PlotRepr.Pcoa, "metric": Metric.Abundance},
            {"Sample", "Classification ID", "PC1", "PC2"},
        ),
        (
            {
                "plot_type": PlotType.Beta,
                "plot_repr": PlotRepr.Distance,
                "metric": Metric.ReadcountWChildren,
            },
            {"Sample", "Sample 1", "Sample 2", "Sample 3"},
        ),
    ],
)
def test_export_chart_data(sample_collection, default_plot_params_payload, params, columns):
    params = PlotParams.model_validate(
        default_plot_params_payload | params | {"export_format": ExportFormat.Csv}
    )

    result = sample_collection.plot(params)
    df = pd.read_csv(io.BytesIO(result.exported_chart_data))

    assert set(df.columns) == columns

    if params.plot_repr == PlotRepr.Distance:
        assert set(df["Sample"]) == columns - {"Sample"}


def test_validate_plot_params_no_functional_profiles(
    sample_collection, default_plot_params_payload
):
    params = PlotParams.model_validate(
        default_plot_params_payload | {"plot_type": PlotType.Functional}
    )

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
        ({"metric": Metric.Readcount}, "Readcount"),
        ({"metric": Metric.ReadcountWChildren}, "Readcount with children"),
        ({"metric": Metric.NormalizedReadcount}, "Normalized readcount"),
        ({"metric": Metric.NormalizedReadcountWChildren}, "Normalized readcount with children"),
        ({"metric": Metric.Abundance}, "Relative abundance"),
        ({"metric": Metric.AbundanceWChildren}, "Relative abundance"),
        ({"facet_by": "Cohort"}, "Cohort"),
    ],
)
def test_get_plot_title(default_plot_params_payload, params, start):
    params = PlotParams.model_validate(
        default_plot_params_payload | params | {"source_name": "My Project"}
    )
    title = get_plot_title(params)
    assert title == f"{start} plot of My Project samples"


CLASSIFICATION_UUIDS = ["0f4ee4ecb3a3412f", "bef0bc57dd7f4c43", "fadd35dd60074a5c"]


@pytest.fixture
def stats_sample_collection() -> SampleCollection:
    samples = []
    for i, (cls_uuid, cohort, site) in enumerate(
        [
            (CLASSIFICATION_UUIDS[0], "A", "X"),
            (CLASSIFICATION_UUIDS[1], "A", "X"),
            (CLASSIFICATION_UUIDS[2], "B", "Y"),
            (CLASSIFICATION_UUIDS[0], "B", "Y"),
        ]
    ):
        samples.append(
            make_sample(
                classification_uuid=cls_uuid,
                filename=f"sample{i}.fastq",
                sample_name=f"Sample {i}",
                cohort=cohort,
                site=site,
            )
        )
    return SampleCollection(samples)


def test_stats_alpha_diversity(stats_sample_collection):
    params = StatsParams(
        source_name="My Source",
        group_by="cohort",
        stats_type=StatsType.AlphaDiversity,
        metric=Metric.ReadcountWChildren,
        rank=Rank.Species,
    )

    result = stats_sample_collection.stats(params)

    assert result.error is None
    assert result.beta_diversity_results is None
    assert isinstance(result.alpha_diversity_results, AlphaDiversityStatsResults)
    assert result.alpha_diversity_results.test == AlphaDiversityStatsTest.Mannwhitneyu
    assert result.alpha_diversity_results.group_sizes == {"A": 2, "B": 2}
    assert result.alpha_diversity_results.sample_size == 4

    assert isinstance(result.plot_results, PlotResults)
    assert result.plot_results.params.plot_type == PlotType.Alpha
    assert result.plot_results.params.plot_repr is None
    assert result.plot_results.params.source_name == params.source_name
    assert result.plot_results.params.metric == params.metric
    assert result.plot_results.params.rank == params.rank
    assert result.plot_results.params.group_by == params.group_by


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_stats_beta_diversity(stats_sample_collection):
    params = StatsParams(
        group_by="cohort",
        stats_type=StatsType.BetaDiversity,
        metric=Metric.ReadcountWChildren,
        rank=Rank.Species,
    )

    result = stats_sample_collection.stats(params)

    assert result.error is None
    assert result.alpha_diversity_results is None
    assert isinstance(result.beta_diversity_results, BetaDiversityStatsResults)
    assert result.beta_diversity_results.test == BetaDiversityStatsTest.Permanova
    assert result.beta_diversity_results.group_sizes == {"A": 2, "B": 2}
    assert result.beta_diversity_results.sample_size == 4
    assert result.beta_diversity_results.num_permutations > 0
    assert result.plot_results is None


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_stats_auto_metric_resolution(stats_sample_collection):
    params = StatsParams(
        group_by="cohort",
        stats_type=StatsType.BetaDiversity,
        rank=Rank.Genus,
    )
    assert params.metric == Metric.Auto

    result = stats_sample_collection.stats(params)

    assert result.error is None
    assert result.params.metric == stats_sample_collection.automatic_metric
    assert result.params.metric != Metric.Auto


def test_stats_secondary_group_by(stats_sample_collection):
    params = StatsParams(
        group_by="cohort",
        secondary_group_by="site",
        stats_type=StatsType.AlphaDiversity,
        metric=Metric.ReadcountWChildren,
        rank=Rank.Species,
    )

    result = stats_sample_collection.stats(params)

    assert result.error is None
    assert result.alpha_diversity_results.group_by_variable == "cohort_site"


def test_stats_alpha_diversity_forwards_shared_params_to_plot(stats_sample_collection):
    params = StatsParams(
        group_by="cohort",
        secondary_group_by="site",
        stats_type=StatsType.AlphaDiversity,
        metric=Metric.ReadcountWChildren,
        rank=Rank.Species,
        alpha_metric=AlphaDiversityMetric.Simpson,
    )

    result = stats_sample_collection.stats(params)

    assert result.plot_results is not None
    plot_params = result.plot_results.params
    assert plot_params.plot_type == PlotType.Alpha
    assert plot_params.plot_repr is None
    assert plot_params.group_by == "cohort"
    assert plot_params.secondary_group_by == "site"
    assert plot_params.metric == params.metric
    assert plot_params.rank == params.rank
    assert plot_params.alpha_metric == params.alpha_metric


def test_stats_filter_by(stats_sample_collection):
    params = StatsParams(
        group_by="cohort",
        stats_type=StatsType.AlphaDiversity,
        metric=Metric.ReadcountWChildren,
        rank=Rank.Species,
        filter_by="cohort",
        filter_value=["A"],
    )

    result = stats_sample_collection.stats(params)

    assert result.error is not None
    assert "at least 2 groups" in result.error


def test_stats_warning_becomes_error(stats_sample_collection):
    params = StatsParams(
        group_by="cohort",
        stats_type=StatsType.AlphaDiversity,
        metric=Metric.AbundanceWChildren,
        rank=Rank.Species,
    )

    result = stats_sample_collection.stats(params)

    assert result.error is not None
    assert result.alpha_diversity_results is None
    assert result.beta_diversity_results is None
    assert result.plot_results is None


def test_stats_invalid_metadata_field(stats_sample_collection):
    params = StatsParams(
        group_by="nonexistent_field",
        stats_type=StatsType.AlphaDiversity,
        metric=Metric.ReadcountWChildren,
        rank=Rank.Species,
    )

    result = stats_sample_collection.stats(params)

    assert result.error is not None
    assert "does not exist" in result.error
    assert result.alpha_diversity_results is None
    assert result.beta_diversity_results is None
    assert result.plot_results is None


@pytest.mark.parametrize(
    "attr",
    ["group_by", "secondary_group_by", "filter_by", "paired_by"],
)
def test_validate_stats_params_invalid_metadata_field(stats_sample_collection, attr):
    kwargs = {"group_by": "cohort", attr: "foo"}
    params = StatsParams(
        stats_type=StatsType.AlphaDiversity,
        metric=Metric.ReadcountWChildren,
        rank=Rank.Species,
        **kwargs,
    )

    with pytest.raises(ValidationError, match="metadata field 'foo' does not exist"):
        stats_sample_collection._validate_stats_params(params)


def test_validate_stats_params_no_classifications():
    samples = [
        Samples(
            {
                "uuid": generate_id(),
                "metadata": {
                    "sample_id": generate_id(),
                    "metadata_id": generate_id(),
                    "classification_id": generate_id(),
                },
                "primary_classification": None,
                "functional_profile": None,
            }
        )
    ]
    collection = SampleCollection(samples)
    params = StatsParams(group_by="cohort", stats_type=StatsType.AlphaDiversity)

    with pytest.raises(ValidationError, match="Classification results are not available"):
        collection._validate_stats_params(params)


def test_stats_to_dict_alpha(stats_sample_collection):
    params = StatsParams(
        group_by="cohort",
        stats_type=StatsType.AlphaDiversity,
        metric=Metric.ReadcountWChildren,
        rank=Rank.Species,
    )

    result = stats_sample_collection.stats(params)
    d = result.to_dict()

    assert d["error"] is None
    assert d["params"]["group_by"] == "cohort"
    assert d["params"]["stats_type"] == "alpha_diversity"
    assert d["beta_diversity_results"] is None
    assert d["alpha_diversity_results"]["test"] == "mannwhitneyu"
    assert isinstance(d["alpha_diversity_results"]["statistic"], float)
    assert isinstance(d["alpha_diversity_results"]["pvalue"], float)
    assert d["alpha_diversity_results"]["group_sizes"] == {"A": 2, "B": 2}
    assert d["alpha_diversity_results"]["paired_by_variable"] is None
    assert "num_permutations" not in d["alpha_diversity_results"]

    assert d["plot_results"] is not None
    assert d["plot_results"]["params"]["plot_type"] == "alpha"
    assert d["plot_results"]["params"]["plot_repr"] is None


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_stats_to_dict_beta(stats_sample_collection):
    params = StatsParams(
        group_by="cohort",
        stats_type=StatsType.BetaDiversity,
        metric=Metric.ReadcountWChildren,
        rank=Rank.Species,
    )

    result = stats_sample_collection.stats(params)
    d = result.to_dict()

    assert d["error"] is None
    assert d["alpha_diversity_results"] is None
    assert d["beta_diversity_results"]["test"] == "permanova"
    assert isinstance(d["beta_diversity_results"]["num_permutations"], int)
    assert "paired_by_variable" not in d["beta_diversity_results"]
    assert d["plot_results"] is None


def test_stats_to_dict_error(stats_sample_collection):
    params = StatsParams(
        group_by="nonexistent_field",
        stats_type=StatsType.AlphaDiversity,
    )

    result = stats_sample_collection.stats(params)
    d = result.to_dict()

    assert d["alpha_diversity_results"] is None
    assert d["beta_diversity_results"] is None
    assert d["plot_results"] is None
    assert d["error"] is not None


def test_stats_to_dict_with_posthoc():
    groups = ["g1", "g2", "g3"]
    adjusted_pvalues = pd.DataFrame(
        [[1.0, 0.02, 0.04], [0.02, 1.0, 0.03], [0.04, 0.03, 1.0]],
        index=groups,
        columns=groups,
    )
    pvalues = pd.DataFrame(
        [[1.0, 0.01, 0.02], [0.01, 1.0, 0.015], [0.02, 0.015, 1.0]],
        index=groups,
        columns=groups,
    )
    statistics = pd.DataFrame(
        [[np.nan, 5.0, 3.0], [5.0, np.nan, 4.0], [3.0, 4.0, np.nan]],
        index=groups,
        columns=groups,
    )

    posthoc = PosthocResults(
        test=PosthocStatsTest.Dunn,
        adjustment_method=AdjustmentMethod.BenjaminiHochberg,
        adjusted_pvalues=adjusted_pvalues,
        pvalues=pvalues,
        statistics=statistics,
    )
    stats_results = AlphaDiversityStatsResults(
        test=AlphaDiversityStatsTest.Kruskal,
        statistic=10.5,
        pvalue=0.005,
        alpha=0.05,
        sample_size=12,
        group_by_variable="cohort",
        group_sizes={"g1": 4, "g2": 4, "g3": 4},
        posthoc=posthoc,
    )
    params = StatsParams(group_by="cohort", stats_type=StatsType.AlphaDiversity)
    result = StatsResults(params=params, alpha_diversity_results=stats_results)

    d = result.to_dict()
    posthoc_dict = d["alpha_diversity_results"]["posthoc"]

    assert posthoc_dict["test"] == "dunn"
    assert posthoc_dict["adjustment_method"] == "benjamini_hochberg"

    assert posthoc_dict["adjusted_pvalues"]["g1"]["g2"] == 0.02
    assert posthoc_dict["adjusted_pvalues"]["g2"]["g1"] == 0.02
    assert posthoc_dict["adjusted_pvalues"]["g1"]["g3"] == 0.04

    assert posthoc_dict["pvalues"]["g1"]["g2"] == 0.01
    assert posthoc_dict["pvalues"]["g2"]["g3"] == 0.015

    assert posthoc_dict["statistics"]["g1"]["g2"] == 5.0
    assert posthoc_dict["statistics"]["g2"]["g3"] == 4.0
    assert np.isnan(posthoc_dict["statistics"]["g1"]["g1"])


def _make_ancombc_main_results():
    index = pd.MultiIndex.from_tuples(
        [
            ("Taxon A", "Intercept"),
            ("Taxon A", "cohort[T.B]"),
            ("Taxon B", "Intercept"),
            ("Taxon B", "cohort[T.B]"),
        ],
        names=["Taxon", "Covariate"],
    )
    return pd.DataFrame(
        {
            "Comparison": [
                "Intercept",
                "cohort: B vs A (reference)",
                "Intercept",
                "cohort: B vs A (reference)",
            ],
            "Log2(FC)": [0.0, 1.5, 0.0, -2.0],
            "SE": [0.1, 0.2, 0.1, 0.3],
            "W": [0.0, 7.5, 0.0, -6.7],
            "pvalue": [1.0, 0.001, 1.0, 0.002],
            "qvalue": [1.0, 0.002, 1.0, 0.004],
            "Signif": [False, True, False, True],
        },
        index=index,
    )


def _make_ancombc_global_results():
    index = pd.Index(["Taxon A", "Taxon B"], name="Taxon")
    return pd.DataFrame(
        {
            "W": [5.0, 3.0],
            "pvalue": [0.01, 0.03],
            "qvalue": [0.02, 0.06],
            "Signif": [True, False],
        },
        index=index,
    )


def test_stats_ancombc(stats_sample_collection):
    params = StatsParams(
        group_by="cohort",
        stats_type=StatsType.Ancombc,
        metric=Metric.ReadcountWChildren,
        rank=Rank.Genus,
    )

    result = stats_sample_collection.stats(params)

    assert result.error is None
    assert result.plot_results is None  # no plot if there's no significant results
    assert result.alpha_diversity_results is None
    assert result.beta_diversity_results is None
    assert isinstance(result.ancombc_results, AncombcResults)
    assert result.ancombc_results.group_sizes == {"A": 2, "B": 2}
    assert result.ancombc_results.sample_size == 4
    assert result.ancombc_results.reference_group == "A"
    assert result.ancombc_results.group_by_variable == "cohort"


def test_stats_ancombc_reference_group(stats_sample_collection):
    params = StatsParams(
        group_by="cohort",
        stats_type=StatsType.Ancombc,
        metric=Metric.ReadcountWChildren,
        rank=Rank.Genus,
        reference_group="B",
    )

    result = stats_sample_collection.stats(params)

    assert result.error is None
    assert result.ancombc_results.reference_group == "B"


def test_stats_ancombc_generates_plot_when_significant(stats_sample_collection):
    main_results = _make_ancombc_main_results()
    ancombc_results = AncombcResults(
        main_results=main_results,
        reference_group="A",
        adjustment_method=AdjustmentMethod.BenjaminiHochberg,
        alpha=0.05,
        sample_size=10,
        group_by_variable="cohort",
        group_sizes={"A": 5, "B": 5},
    )
    params = StatsParams(
        group_by="cohort",
        stats_type=StatsType.Ancombc,
        metric=Metric.ReadcountWChildren,
        rank=Rank.Species,
    )
    mock_stats_results = StatsResults(params=params, ancombc_results=ancombc_results)

    with mock.patch.object(
        type(stats_sample_collection), "_stats", return_value=mock_stats_results
    ):
        result = stats_sample_collection.stats(params)

    assert result.error is None
    assert isinstance(result.plot_results, PlotResults)
    assert result.plot_results.error is None
    assert isinstance(result.plot_results.chart, dict)


def test_stats_to_dict_ancombc():
    main_results = _make_ancombc_main_results()
    global_results = _make_ancombc_global_results()
    ancombc_results = AncombcResults(
        main_results=main_results,
        global_results=global_results,
        reference_group="A",
        adjustment_method=AdjustmentMethod.BenjaminiHochberg,
        alpha=0.05,
        sample_size=10,
        group_by_variable="cohort",
        group_sizes={"A": 5, "B": 5},
    )
    params = StatsParams(group_by="cohort", stats_type=StatsType.Ancombc)
    result = StatsResults(params=params, ancombc_results=ancombc_results)

    d = result.to_dict()

    assert d["error"] is None
    assert d["alpha_diversity_results"] is None
    assert d["beta_diversity_results"] is None
    assert d["params"]["stats_type"] == "ancombc"

    ancombc = d["ancombc_results"]
    assert ancombc["reference_group"] == "A"
    assert ancombc["adjustment_method"] == "benjamini_hochberg"
    assert ancombc["alpha"] == 0.05
    assert ancombc["sample_size"] == 10
    assert ancombc["group_by_variable"] == "cohort"
    assert ancombc["group_sizes"] == {"A": 5, "B": 5}

    records = ancombc["main_results"]
    assert isinstance(records, list)
    assert len(records) == 4
    assert all(
        {"Taxon", "Covariate", "Comparison", "Log2(FC)", "SE", "W", "pvalue", "qvalue", "Signif"}
        == set(r.keys())
        for r in records
    )
    sig_record = next(
        r for r in records if r["Taxon"] == "Taxon A" and r["Covariate"] == "cohort[T.B]"
    )
    assert sig_record["Log2(FC)"] == 1.5
    assert sig_record["Signif"] is True

    global_records = ancombc["global_results"]
    assert isinstance(global_records, list)
    assert len(global_records) == 2
    assert all(
        {"Taxon", "W", "pvalue", "qvalue", "Signif"} == set(r.keys()) for r in global_records
    )
    taxon_a = next(r for r in global_records if r["Taxon"] == "Taxon A")
    assert taxon_a["W"] == 5.0
    assert taxon_a["Signif"] is True
