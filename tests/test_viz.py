import math
import mock
import pytest
from datetime import datetime

pytest.importorskip("numpy")  # noqa
pytest.importorskip("pandas")  # noqa

import altair as alt
import numpy as np
import pandas as pd

from onecodex.lib.enums import Metric, AbundanceMetric, Link, Rank, BetaDiversityMetric, Linkage
from onecodex.exceptions import OneCodexException, PlottingException
from onecodex.models.collection import SampleCollection
from onecodex.utils import has_missing_values
from onecodex.viz._primitives import (
    interleave_palette,
    get_classification_url,
    get_ncbi_taxonomy_browser_url,
)
from onecodex.exceptions import PlottingWarning


def assert_expected_legend_title(chart, title):
    assert chart.to_dict()["encoding"]["color"]["legend"]["title"] == title


def assert_expected_x_axis_title(chart, title):
    assert chart.to_dict()["encoding"]["x"]["axis"]["title"] == title


def assert_expected_y_axis_title(chart, title):
    assert chart.to_dict()["encoding"]["y"]["axis"]["title"] == title


def assert_expected_color_scale_domain_len(chart, length):
    assert len(chart.to_dict()["encoding"]["color"]["scale"]["domain"]) == length


def assert_expected_color_scale_range_len(chart, length):
    assert len(chart.to_dict()["encoding"]["color"]["scale"]["range"]) == length


def assert_expected_color_field(chart, field):
    assert chart.to_dict()["encoding"]["color"]["field"] == field


def assert_expected_x_offset_field(chart, field):
    assert chart.to_dict()["encoding"]["xOffset"]["field"] == field


def assert_expected_column_field(chart, field):
    assert chart.to_dict()["encoding"]["column"]["field"] == field


def test_altair_ocx_theme(ocx, api_data):
    assert alt.theme.active == "onecodex"


def test_altair_ocx_renderer(ocx, api_data):
    assert alt.renderers.active == "onecodex"


def test_plot_metadata(samples):
    chart = samples.plot_metadata(
        vaxis="simpson", title="my title", xlabel="my xlabel", ylabel="my ylabel", return_chart=True
    )
    assert chart.data["simpson"].tolist() == [
        0.9232922257199748,
        0.8930761430647977,
        0.7865654458730156,
    ]
    assert chart.mark == "circle"
    assert chart.title == "my title"
    assert_expected_x_axis_title(chart, "my xlabel")
    assert_expected_y_axis_title(chart, "my ylabel")

    # try time, boolean, and numerical types for x-axis
    chart = samples.plot_metadata(haxis="date_sequenced", vaxis="observed_taxa", return_chart=True)
    assert chart.encoding.x.shorthand == "date_sequenced"

    chart = samples.plot_metadata(haxis="starred", vaxis="observed_taxa", return_chart=True)
    assert chart.encoding.x.shorthand == "starred"

    chart = samples.plot_metadata(haxis="totalige", vaxis="observed_taxa", return_chart=True)
    assert chart.mark == "circle"
    assert chart.encoding.x.shorthand == "totalige"

    # taxid and taxon on vertical axis
    chart = samples.plot_metadata(vaxis=853, plot_type="boxplot", return_chart=True)
    assert chart.mark.type == "boxplot"

    chart = samples.plot_metadata(
        vaxis="Faecalibacterium prausnitzii", plot_type="scatter", return_chart=True
    )
    assert chart.mark == "circle"


def test_plot_metadata_facet_by_scatter(samples):
    chart = samples.plot_metadata(vaxis="shannon", facet_by="geo_loc_name", return_chart=True)

    assert chart.mark == "circle"
    assert_expected_x_axis_title(chart, "")


def test_plot_metadata_facet_by_boxplot(samples):
    chart = samples.plot_metadata(
        vaxis="shannon", haxis="starred", facet_by="geo_loc_name", return_chart=True
    )

    assert chart.mark.type == "boxplot"
    assert_expected_x_axis_title(chart, "")


def test_plot_metadata_secondary_haxis_scatter(samples):
    chart = samples.plot_metadata(
        vaxis="shannon", haxis="geo_loc_name", secondary_haxis="Label", return_chart=True
    )

    assert chart.mark == "circle"
    assert_expected_x_axis_title(chart, "geo_loc_name")
    assert_expected_color_field(chart, "Label")
    assert_expected_x_offset_field(chart, "Label")


def test_plot_metadata_secondary_haxis_boxplot(samples):
    chart = samples.plot_metadata(
        vaxis="shannon",
        haxis="geo_loc_name",
        secondary_haxis="starred",
        plot_type="boxplot",
        return_chart=True,
    )

    assert chart.mark.type == "boxplot"
    assert_expected_x_axis_title(chart, "geo_loc_name")
    assert_expected_color_field(chart, "starred")
    assert_expected_x_offset_field(chart, "starred")


@pytest.mark.parametrize("plot_type", ["scatter", "boxplot"])
def test_plot_metadata_secondary_haxis_and_facet_by(samples, plot_type):
    chart = samples.plot_metadata(
        vaxis="shannon",
        haxis="geo_loc_name",
        secondary_haxis="Label",
        facet_by="starred",
        plot_type=plot_type,
        return_chart=True,
    )

    if plot_type == "scatter":
        assert chart.mark == "circle"
    else:
        assert chart.mark.type == "boxplot"

    assert_expected_x_axis_title(chart, "")
    assert_expected_color_field(chart, "Label")
    assert_expected_x_offset_field(chart, "Label")
    assert_expected_column_field(chart, "starred")


def test_plot_metadata_alpha_diversity_with_nans(samples):
    mock_alpha_div_values = samples.alpha_diversity(metric="shannon")
    mock_alpha_div_values.iat[1, 0] = np.nan
    mock_alpha_div_values.iat[2, 0] = np.nan

    with mock.patch.object(SampleCollection, "alpha_diversity", return_value=mock_alpha_div_values):
        chart = samples.plot_metadata(vaxis="shannon", return_chart=True)

    shannon_result = chart.data["shannon"].tolist()
    assert len(shannon_result) == 1
    assert math.isclose(shannon_result[0], 3.6127134403472887)
    assert chart.mark == "circle"


@pytest.mark.parametrize("coerce_haxis_dates", [True, False])
def test_plot_metadata_haxis_date_coercion(samples, coerce_haxis_dates):
    samples.metadata["date_collected"] = datetime.now().isoformat()
    assert samples.metadata["date_collected"].dtype == object

    chart = samples.plot_metadata(
        vaxis="shannon",
        haxis="date_collected",
        return_chart=True,
        coerce_haxis_dates=coerce_haxis_dates,
    )

    assert chart.mark.type == "boxplot"
    assert chart.encoding.x.shorthand == "date_collected"
    # previous column type coercion would raise a json serialization TypeError
    assert isinstance(chart.to_dict(), dict)

    if coerce_haxis_dates:
        assert pd.api.types.is_datetime64_any_dtype(chart.data["date_collected"])
    else:
        assert chart.data["date_collected"].dtype == object


def test_plot_metadata_exceptions(samples):
    # expect error if rank is None, since that could lead to weird results
    with pytest.raises(OneCodexException) as e:
        samples.plot_metadata(vaxis="simpson", rank=None)
    assert "specify a rank" in str(e.value)

    # vert axis does not exist
    with pytest.raises(OneCodexException):
        samples.plot_metadata(vaxis="does_not_exist")

    # must have at least one sample
    with pytest.raises(PlottingException) as e:
        samples[:0].plot_metadata()
    assert "too few samples" in str(e.value)

    with pytest.raises(OneCodexException, match="`haxis` and `secondary_haxis` cannot be the same"):
        samples.plot_metadata(vaxis="shannon", haxis="wheat", secondary_haxis="wheat")


def test_plot_metadata_all_samples_are_nan(samples):
    samples._results[:] = np.nan
    assert len(samples._classification_ids_without_abundances) == 3

    # Raises for alpha diversity
    with pytest.raises(PlottingException, match="Abundances are not calculated for any"):
        samples.plot_metadata(vaxis="shannon")

    # Does not raise for numeric metadata column
    chart = samples.plot_metadata(vaxis="totalige", return_chart=True)
    assert chart.mark == "circle"


def test_plot_metadata_filters_nan_samples(samples):
    samples._results.iloc[0, :] = np.nan
    samples._results.iloc[1, :] = np.nan
    assert len(samples._classification_ids_without_abundances) == 2

    with pytest.warns(PlottingWarning, match=r"2 sample\(s\) have no abundances calculated"):
        chart = samples.plot_metadata(vaxis="shannon", return_chart=True)

    assert len(chart["data"]) == 1


def test_plot_metadata_group_with_single_value(samples):
    chart = samples.plot_metadata(
        plot_type="boxplot",
        haxis=("library_type", "external_sample_id"),
        vaxis="observed_taxa",
        return_chart=True,
    )
    assert chart.mark.type == "boxplot"

    chart = samples.plot_metadata(
        plot_type="boxplot",
        haxis="library_type",
        vaxis="observed_taxa",
        facet_by="external_sample_id",
        return_chart=True,
    )
    assert chart.mark.type == "boxplot"


def test_plot_pca(samples):
    samples._collate_results(metric="readcount_w_children")

    chart = samples.plot_pca(
        title="my title",
        xlabel="my xlabel",
        ylabel="my ylabel",
        color="geo_loc_name",
        size="totalige",
        rank="genus",
        org_vectors=3,
        return_chart=True,
        tooltip=["totalige", "vegetables", "Prevotella", 816],
    )

    assert chart.title == "my title"

    # one for main plot, one for eigenvectors
    assert len(chart.layer) == 2

    mainplot = chart.layer[0]
    assert mainplot.data["Bacteroides (816)"].round(6).tolist() == [0.356439, 0.213307, 0.787949]
    assert mainplot.data["Prevotella (838)"].round(6).tolist() == [0.000221, 0.001872, 6.7e-05]
    assert mainplot.data["totalige"].tolist() == [62.9, 91.5, 112.0]
    assert mainplot.encoding.color.shorthand == "geo_loc_name"
    assert mainplot.encoding.size.shorthand == "totalige"
    assert mainplot.encoding.x.shorthand == "PC1"
    assert_expected_x_axis_title(mainplot, "my xlabel")
    assert mainplot.encoding.y.shorthand == "PC2"
    assert_expected_y_axis_title(mainplot, "my ylabel")
    assert sorted([x.shorthand for x in mainplot.encoding.tooltip]) == [
        "Bacteroides (816)",
        "Label",
        "Prevotella (838)",
        "geo_loc_name",
        "totalige",
        "vegetables",
    ]

    vectors = chart.layer[1]
    assert vectors.data["x"].sum().round(6) == 0.145172
    assert vectors.data["y"].sum().round(6) == 0.039944
    assert vectors.data["o"].tolist() == [0, 1, 0, 1, 0, 1]


@pytest.mark.parametrize("match_taxonomy", [True, False])
def test_plot_match_taxonomy(samples, match_taxonomy):
    """
    Test that arguments that determine plotting parameters based on
    metadata/taxonomy can be configured to exclude taxonomy
    """
    chart = samples.plot_pca(
        tooltip="Bacteroides", match_taxonomy=match_taxonomy, return_chart=True
    )
    assert ("Bacteroides (816)" in chart.data.columns) == match_taxonomy

    chart = samples.plot_heatmap(
        tooltip="Bacteroides", match_taxonomy=match_taxonomy, return_chart=True
    )
    assert ("Bacteroides (816)" in chart.data.columns) == match_taxonomy

    chart = samples.plot_bargraph(
        tooltip="Bacteroides", match_taxonomy=match_taxonomy, return_chart=True
    )
    assert ("Bacteroides (816)" in chart.data.columns) == match_taxonomy

    # returns an HConcatChart
    chart = samples.plot_distance(
        tooltip="Bacteroides", match_taxonomy=match_taxonomy, return_chart=True
    )
    # `plot_distance` returns an HConcatChart, so a little digging is required
    assert (
        "2) Bacteroides (816)" in list(chart._kwds["hconcat"][1].data.columns)
    ) == match_taxonomy

    # plot_metadata does not have a tooltip argument, so test with haxis instead
    chart = samples.plot_metadata(
        haxis="Bacteroides", match_taxonomy=match_taxonomy, return_chart=True
    )
    assert ("Bacteroides (816)" in chart.data.columns) == match_taxonomy

    chart = samples.plot_mds(
        tooltip="Bacteroides", match_taxonomy=match_taxonomy, return_chart=True
    )
    assert ("Bacteroides (816)" in chart.data.columns) == match_taxonomy

    chart = samples.plot_pcoa(
        tooltip="Bacteroides", match_taxonomy=match_taxonomy, return_chart=True
    )
    assert ("Bacteroides (816)" in chart.data.columns) == match_taxonomy


def test_plot_pca_color_by_bool_field(samples):
    assert samples.metadata["wheat"].dtype == bool

    chart = samples.plot_pca(color="wheat", return_chart=True)

    assert chart.encoding.color.shorthand == "wheat"
    assert_expected_legend_title(chart, "wheat")
    assert_expected_color_scale_domain_len(chart, 3)
    assert_expected_color_scale_range_len(chart, 2)


def test_plot_pca_color_by_field_with_nans(samples):
    assert has_missing_values(samples.metadata["name"])

    chart = samples.plot_pca(color="name", return_chart=True)

    assert chart.encoding.color.shorthand == "name"
    assert_expected_legend_title(chart, "name")
    assert_expected_color_scale_domain_len(chart, 3)
    assert_expected_color_scale_range_len(chart, 1)


def test_plot_pca_missing_abundances(ocx, api_data, samples):
    sample1 = ocx.Samples.get("cc18208d98ad48b3")
    sample2 = ocx.Samples.get("5445740666134eee")
    samples.extend([sample1, sample2])

    samples._results.loc[samples[1].primary_classification.id] = np.nan
    samples._results.loc[samples[3].primary_classification.id] = np.nan
    assert len(samples._classification_ids_without_abundances) == 2

    with pytest.warns(PlottingWarning, match=r"2 sample\(s\) have no abundances calculated"):
        samples.plot_pca(return_chart=True)


def test_plot_pca_exceptions(samples):
    # expect error if rank is None, since that could lead to weird results
    with pytest.raises(OneCodexException) as e:
        samples.plot_pca(rank=None)
    assert "specify a rank" in str(e.value)

    # must have at least three samples
    with pytest.raises(PlottingException) as e:
        samples[:2].plot_pca()
    assert "too few samples" in str(e.value)

    # samples must have at least two taxa at this rank
    with pytest.raises(PlottingException) as e:
        samples.to_df(top_n=1).ocx.plot_pca()
    assert "too few taxa" in str(e.value)


def test_plot_heatmap(samples):
    chart = samples.plot_heatmap(
        top_n=10, title="my title", xlabel="my xlabel", ylabel="my ylabel", return_chart=True
    )
    assert chart.mark == "rect"
    assert chart.title == "my title"
    assert_expected_x_axis_title(chart, "my xlabel")
    assert_expected_y_axis_title(chart, "my ylabel")
    assert len(chart.data["tax_id"].unique()) == 10
    assert chart.data["Relative Abundance"].sum().round(6) == 1.813541

    chart = samples.plot_heatmap(threshold=0.1, haxis="eggs", return_chart=True)
    assert chart.mark == "rect"
    assert all(chart.data.groupby("tax_id").max()["Relative Abundance"] > 0.1)

    chart = samples.plot_heatmap(top_n=10, threshold=0.01, return_chart=True)
    assert len(chart.data["tax_id"].unique()) == 10
    assert all(chart.data.groupby("tax_id").max()["Relative Abundance"] > 0.01)


def test_plot_heatmap_match_taxonomy(samples):
    chart = samples.plot_heatmap(
        top_n=10, tooltip="Bacteroides", match_taxonomy=True, return_chart=True
    )
    assert "Bacteroides (816)" in chart.data.columns

    chart = samples.plot_heatmap(
        top_n=10, tooltip="Bacteroides", match_taxonomy=False, return_chart=True
    )
    assert "Bacteroides (816)" not in chart.data.columns


def test_plot_heatmap_with_missing_haxis_sample(samples):
    samples[2].metadata.custom.pop("eggs")

    # Does not raise exception if a sample is missing `haxis` value in custom metadata
    chart = samples.plot_heatmap(top_n=10, threshold=0.1, haxis="eggs", return_chart=True)

    assert "N/A" in chart.data["eggs"].unique()


@pytest.mark.parametrize(
    "is_onecodex_accessor",
    [
        (False),
        (True),
    ],
)
def test_plot_heatmap_plots_samples_without_abundances_with_nans(samples, is_onecodex_accessor):
    # Set abundances for first sample to be all NaN
    all_nan_sample = samples[0]
    samples._results.loc[all_nan_sample.primary_classification.id] = np.nan
    assert len(samples._classification_ids_without_abundances) == 1

    all_nan_classification_id = samples._classification_ids_without_abundances[0]
    if is_onecodex_accessor:
        chart = samples.to_df().ocx.plot_heatmap(
            title="my title", xlabel="my xlabel", ylabel="my ylabel", return_chart=True
        )
    else:
        chart = samples.plot_heatmap(
            top_n=10, title="my title", xlabel="my xlabel", ylabel="my ylabel", return_chart=True
        )
    chart_data = chart["data"]

    # Samples with no abundances calculated for them should have all-NaN rows in chart data df
    # Other classification_id rows should not have any `NaN`s
    assert (
        chart_data[(chart_data == all_nan_classification_id).any(axis=1)]["Relative Abundance"]
        .isnull()
        .all()
        == np.True_
    )
    assert (
        chart_data[(chart_data == samples[1].primary_classification.id).any(axis=1)][
            "Relative Abundance"
        ]
        .isnull()
        .any()
        == np.False_
    )

    # Sample with all NaNs should be at the end
    assert chart.to_dict()["encoding"]["x"]["sort"][-1] == all_nan_sample.filename


def test_plot_heatmap_with_haxis_two_cluster_groups(ocx, api_data, samples):
    sample1 = ocx.Samples.get("cc18208d98ad48b3")
    sample2 = ocx.Samples.get("5445740666134eee")
    samples.extend([sample1, sample2])

    # Set abundances for one sample to be all NaN
    # (It should not be in the clustering group)
    all_nan_sample = samples[1]
    samples._results.loc[all_nan_sample.primary_classification.id] = np.nan
    assert len(samples._classification_ids_without_abundances) == 1

    # Make sure we are actually calling _cluster_by_sample
    with mock.patch("onecodex.viz._distance.VizDistanceMixin._cluster_by_sample") as cluster_fn:
        # The _cluster_by_sample method returns other things as well, but `ids_in_order` is the
        # only thing needed for the rest of `plot_heatmap`
        cluster_fn.return_value = {
            "ids_in_order": [
                "4eed25415f6945dd",
                "6579e99943f84ad2",
                "e0422602de41479f",
                "1198c7ce565643f9",
            ]
        }
        samples.plot_heatmap(rank="genus", top_n=15, haxis="wheat", return_chart=True)
        assert cluster_fn.call_count == 1

    # Does not raise exception
    samples.plot_heatmap(rank="genus", top_n=15, haxis="wheat", return_chart=True)


def test_plot_heatmap_all_samples_are_nan(ocx, api_data, samples):
    samples._results[:] = np.nan
    assert len(samples._classification_ids_without_abundances) == 3

    with pytest.raises(PlottingException, match="Abundances are not calculated for any"):
        samples.plot_heatmap(top_n=10, return_chart=True)


def test_plot_bargraph_no_samples_have_abundances(ocx, api_data, samples):
    samples._results[:] = np.nan
    assert len(samples._classification_ids_without_abundances) == 3

    chart = samples.plot_bargraph(return_chart=True)

    assert list(chart.data["Relative Abundance"].values) == [0.0, 0.0, 0.0]


def test_plot_distance_cluster_by_sample_excludes_classifications_without_abundances(
    ocx, api_data, samples
):
    sample1 = ocx.Samples.get("cc18208d98ad48b3")
    sample2 = ocx.Samples.get("5445740666134eee")
    samples.extend([sample1, sample2])

    # Set abundances for one sample to be all NaN
    all_nan_sample = samples[1]
    samples._results.loc[all_nan_sample.primary_classification.id] = np.nan
    assert len(samples._classification_ids_without_abundances) == 1

    with (
        mock.patch("onecodex.viz._distance.VizDistanceMixin._cluster_by_sample") as cluster_fn,
        mock.patch("onecodex.viz.dendrogram"),
        mock.patch("altair.hconcat"),
    ):
        with pytest.warns(PlottingWarning):
            samples.plot_distance()
            # the cluster function should be called with
            # `exclude_classifications_without_abundances=True` and the classification ID of the
            # `all_nan_sample`
            cluster_fn.assert_called_with(
                rank=Rank.Auto,
                metric=BetaDiversityMetric.BrayCurtis,
                linkage=Linkage.Average,
                classification_ids_without_abundances=[all_nan_sample.primary_classification.id],
                exclude_classifications_without_abundances=True,
            )


def test_plot_distance_excludes_classifications_without_abundances(ocx, api_data, samples):
    sample1 = ocx.Samples.get("cc18208d98ad48b3")
    sample2 = ocx.Samples.get("5445740666134eee")
    samples.extend([sample1, sample2])

    # Set abundances for one sample to be all NaN
    all_nan_sample = samples[1]
    samples._results.loc[all_nan_sample.primary_classification.id] = np.nan
    assert len(samples._classification_ids_without_abundances) == 1

    with mock.patch("altair.hconcat") as mock_hconcat:
        with pytest.warns(PlottingWarning):
            samples.plot_distance()
            assert (
                all_nan_sample.primary_classification.id
                not in mock_hconcat.call_args[0][1].data["1) Label"].values
            )


def test_plot_distance_min_with_abundances(ocx, api_data):
    import numpy as np

    sample1 = ocx.Samples.get("cc18208d98ad48b3")
    sample2 = ocx.Samples.get("5445740666134eee")
    samples = SampleCollection([sample1, sample2])

    # Set abundances for both samples to be all NaN
    samples._results.loc[sample1.primary_classification.id] = np.nan
    samples._results.loc[sample2.primary_classification.id] = np.nan
    assert len(samples._classification_ids_without_abundances) == 2

    # We should raise a PlottingException if we don't have >= 2 samples with abundances calculated
    with pytest.raises(PlottingException) as e:
        samples.plot_distance()
        assert "There are too few samples for distance matrix plots" in str(e.value)


def test_plot_heatmap_exceptions(samples):
    # expect error if rank is None, since that could lead to weird results
    with pytest.raises(OneCodexException) as e:
        samples.plot_heatmap(rank=None)
    assert "specify a rank" in str(e.value)

    # must have at least two samples
    with pytest.raises(PlottingException) as e:
        samples[:1].plot_heatmap()
    assert "too few samples" in str(e.value)

    # must have at least two taxa
    with pytest.raises(PlottingException) as e:
        samples.plot_heatmap(top_n=1)
    assert "too few taxa" in str(e.value)

    # must specify at least threshold or top_n
    with pytest.raises(OneCodexException) as e:
        samples.plot_heatmap(top_n=None, threshold=None)
    assert "specify at least one of" in str(e.value)


def test_plot_distance(samples):
    chart = samples.plot_distance(
        metric="weighted_unifrac",
        xlabel="my xlabel",
        ylabel="my ylabel",
        title="my title",
        tooltip="vegetables",
        return_chart=True,
    )
    assert len(chart.hconcat) == 2

    mainplot = chart.hconcat[1]
    assert mainplot.mark == "rect"
    assert_expected_x_axis_title(mainplot, "my xlabel")
    assert_expected_y_axis_title(mainplot, "my ylabel")
    assert sorted([x.shorthand for x in mainplot.encoding.tooltip]) == [
        "1) Label",
        "1) vegetables",
        "2) Label",
        "2) vegetables",
        "Distance:Q",
    ]
    assert mainplot.data["Distance"].sum() == 3.022956


def test_plot_distance_exceptions(samples):
    # expect error if rank is None, since that could lead to weird results
    with pytest.raises(OneCodexException) as e:
        samples.plot_distance(rank=None)
    assert "specify a rank" in str(e.value)

    # only some metrics allowed
    with pytest.raises(OneCodexException) as e:
        samples.plot_distance(metric="simpson")
    assert "must be one of" in str(e.value)

    # need at least two samples
    with pytest.raises(PlottingException) as e:
        samples[:1].plot_distance(
            metric="jaccard", xlabel="my xlabel", ylabel="my ylabel", title="my title"
        )
    assert "too few samples" in str(e.value)


@pytest.mark.parametrize(
    "metric,dissimilarity_metric,smacofs",
    [
        # Some results differ between sklearn 1.6 and 1.7
        ("abundance_w_children", "weighted_unifrac", {0.7595}),
        ("abundance_w_children", "unweighted_unifrac", {0.1734}),
        ("abundance_w_children", "braycurtis", {0.0143, 0.1284}),
        ("readcount_w_children", "weighted_unifrac", {0.4956, 0.0918}),
        ("readcount_w_children", "unweighted_unifrac", {0.3579}),
        ("readcount_w_children", "braycurtis", {0.1735, 0.0798}),
    ],
)
def test_plot_mds(samples, metric, dissimilarity_metric, smacofs):
    samples._collate_results(metric=metric)

    chart = samples.plot_mds(
        rank="species",
        method="pcoa",
        metric=dissimilarity_metric,
        xlabel="my xlabel",
        ylabel="my ylabel",
        title="my title",
        return_chart=True,
    )
    assert chart.mark.type == "circle"
    assert chart.title == "my title"
    assert chart.encoding.x.shorthand == "PC1"
    assert_expected_x_axis_title(chart, "my xlabel")
    assert chart.encoding.y.shorthand == "PC2"
    assert_expected_y_axis_title(chart, "my ylabel")
    assert (chart.data["PC1"] * chart.data["PC2"]).sum().round(6) == 0.0

    chart = samples.plot_mds(
        method="smacof", rank="species", metric=dissimilarity_metric, return_chart=True
    )
    assert (chart.data["MDS1"] * chart.data["MDS2"]).sum().round(4) in smacofs


def test_plot_mds_color_by_bool_field(samples):
    assert samples.metadata["wheat"].dtype == bool

    chart = samples.plot_mds(color="wheat", return_chart=True)

    assert chart.encoding.color.shorthand == "wheat"
    assert_expected_legend_title(chart, "wheat")
    assert_expected_color_scale_domain_len(chart, 3)
    assert_expected_color_scale_range_len(chart, 2)


def test_plot_mds_color_by_field_with_nans(samples):
    assert has_missing_values(samples.metadata["name"])

    chart = samples.plot_mds(color="name", return_chart=True)

    assert chart.encoding.color.shorthand == "name"
    assert_expected_legend_title(chart, "name")
    assert_expected_color_scale_domain_len(chart, 3)
    assert_expected_color_scale_range_len(chart, 1)


def test_plot_mds_missing_abundances(ocx, api_data, samples):
    sample1 = ocx.Samples.get("cc18208d98ad48b3")
    sample2 = ocx.Samples.get("5445740666134eee")
    samples.extend([sample1, sample2])

    samples._results.loc[samples[1].primary_classification.id] = np.nan
    samples._results.loc[samples[3].primary_classification.id] = np.nan
    assert len(samples._classification_ids_without_abundances) == 2

    with pytest.warns(PlottingWarning, match=r"2 sample\(s\) have no abundances calculated"):
        samples.plot_mds(return_chart=True)


def test_plot_pcoa(samples):
    chart = samples.plot_pcoa(
        metric="weighted_unifrac",
        xlabel="my xlabel",
        ylabel="my ylabel",
        title="my title",
        return_chart=True,
    )

    assert (chart.data["PC1"] * chart.data["PC2"]).sum().round(6) == 0.0


def test_plot_pcoa_exceptions(samples):
    # need at least 3 samples
    with pytest.raises(PlottingException) as e:
        samples[:2].plot_pcoa()
    assert "too few samples" in str(e.value)


def test_plot_mds_exceptions(samples):
    # expect error if rank is None, since that could lead to weird results
    with pytest.raises(OneCodexException) as e:
        samples.plot_mds(rank=None)
    assert "specify a rank" in str(e.value)

    # only some metrics allowed
    with pytest.raises(OneCodexException) as e:
        samples.plot_mds(metric="simpson")
    assert "must be one of" in str(e.value)

    # need at least 3 samples
    with pytest.raises(PlottingException) as e:
        samples[:2].plot_mds(
            metric="jaccard", xlabel="my xlabel", ylabel="my ylabel", title="my title"
        )
    assert "too few samples" in str(e.value)


@pytest.mark.parametrize(
    "kwargs,metric,msg",
    [
        ({"rank": None}, None, "specify a rank"),
        ({"threshold": None, "top_n": None}, None, "threshold, top_n"),
        ({"normalize": False}, Metric.AbundanceWChildren, "already been normalized"),
        ({"group_by": "foo", "tooltip": "bar"}, None, "not supported with `group_by`"),
        ({"group_by": "foo", "haxis": "bar"}, None, "not supported with `group_by`"),
        ({"group_by": "foo", "label": "bar"}, None, "not supported with `group_by`"),
        ({"group_by": "foo", "sort_x": []}, None, "not supported with `group_by`"),
        ({"group_by": "foo"}, None, "field foo not found"),
    ],
)
def test_plot_bargraph_errors(samples, kwargs, metric, msg):
    if metric:
        samples._collate_results(metric=metric)

    with pytest.raises(OneCodexException, match=msg):
        samples.plot_bargraph(**kwargs)


def test_plot_bargraph_too_few_samples(samples):
    # need at least 1 sample
    with pytest.raises(PlottingException, match="too few samples"):
        samples[:0].plot_bargraph()


@pytest.mark.parametrize("metric", [Metric.Readcount, Metric.ReadcountWChildren])
def test_plot_bargraph_group_by_counts_already_normalized(samples, metric):
    samples._collate_results(metric=metric)

    with pytest.raises(OneCodexException, match="`group_by`.*readcounts.*already been normalized"):
        with mock.patch.object(samples, "_guess_normalized", return_value=True):
            samples.plot_bargraph(group_by="barley")


@pytest.mark.parametrize(
    "metric,rank,label",
    [
        ("abundance_w_children", "genus", "Relative Abundance"),
        ("abundance_w_children", "species", "Relative Abundance"),
        ("readcount_w_children", "genus", "Reads"),
        ("readcount_w_children", "species", "Reads"),
    ],
)
def test_plot_bargraph_chart_result(samples, metric, rank, label):
    samples._collate_results(metric=metric)
    chart = samples.plot_bargraph(
        rank=rank,
        return_chart=True,
        title="Glorious Bargraph",
        xlabel="Exemplary Samples",
        ylabel="Glorious Abundances",
        sort_x=sorted,
        width=200,
        height=200,
    )

    assert chart.title == "Glorious Bargraph"
    assert chart.encoding.x.shorthand == "Label"
    assert_expected_x_axis_title(chart, "Exemplary Samples")
    assert chart.encoding.y.shorthand == label
    assert_expected_y_axis_title(chart, "Glorious Abundances")
    assert_expected_legend_title(chart, label)
    assert chart.encoding.href.shorthand == "url:N"


@pytest.mark.parametrize(
    "metric,kwargs",
    [
        (Metric.Readcount, {}),
        (Metric.ReadcountWChildren, {}),
        (Metric.ReadcountWChildren, {"threshold": 0.1, "top_n": 15, "normalize": True}),
        (Metric.ReadcountWChildren, {"threshold": 0.1, "top_n": None, "normalize": True}),
        (Metric.ReadcountWChildren, {"threshold": None, "top_n": 12, "normalize": True}),
        (Metric.ReadcountWChildren, {"threshold": 2, "top_n": 12, "normalize": False}),
        (Metric.Abundance, {}),
        (Metric.AbundanceWChildren, {}),
    ],
)
def test_plot_bargraph_with_group_by(samples, metric, kwargs):
    samples._collate_results(metric=metric)
    chart = samples.plot_bargraph(group_by="barley", return_chart=True, **kwargs)

    assert chart.encoding.x.shorthand == "barley"
    assert_expected_x_axis_title(chart, "barley")

    metric_label = (
        "Mean Relative Abundance:Q" if AbundanceMetric.has_value(metric) else "Mean Reads:Q"
    )
    assert [x.shorthand for x in chart.encoding.tooltip] == ["barley", "tax_name", metric_label]
    assert not hasattr(chart.encoding.href, "shorthand")


def test_plot_bargraph_include_other_with_empty_samples(samples):
    classification_id = samples._results.index[1]
    samples._results.loc[classification_id] = 0.0

    df = samples.plot_bargraph(rank=Rank.Genus, include_other=True, return_chart=True).data
    total_abundance = df[df["classification_id"] == classification_id]["Relative Abundance"].sum()

    assert total_abundance == 0.0


@pytest.mark.parametrize(
    "legend,expected_title",
    [
        ("auto", "Reads"),
        ("my plot title", "my plot title"),
        (alt.Legend(title="a different title"), "a different title"),
    ],
)
def test_plot_bargraph_legend(samples, legend, expected_title):
    samples._collate_results(metric="readcount_w_children")
    chart = samples.plot_bargraph(return_chart=True, legend=legend)

    assert_expected_legend_title(chart, expected_title)


def test_plot_bargraph_legend_error(samples):
    with pytest.raises(TypeError, match="legend.*list"):
        samples.plot_bargraph(legend=[4, 2])


@pytest.mark.parametrize(
    "plot_type",
    ["bargraph", "heatmap"],
)
@pytest.mark.parametrize(
    "link",
    [Link.Ocx, Link.Ncbi],
)
def test_plot_bargraph_and_heatmap_link(samples, plot_type, link):
    if plot_type == "bargraph":
        chart = samples.plot_bargraph(
            return_chart=True,
            link=link,
        )
    elif plot_type == "heatmap":
        chart = samples.plot_heatmap(
            return_chart=True,
            link=link,
        )
    else:
        raise NotImplementedError()

    assert chart.encoding.href.shorthand == "url:N"

    urls = chart.data["url"]
    for url in urls:
        if link == Link.Ocx:
            assert "/classification/" in url
        elif link == Link.Ncbi:
            if url != "":
                assert url.startswith("https://www.ncbi")
                assert "id=" in url
        else:
            raise NotImplementedError()


def test_plot_bargraph_when_passing_tax_id_as_group_by(samples):
    samples._collate_results()
    for i, sample in enumerate(samples):
        sample.metadata.custom["tax_id"] = f"TAX_{i}"
    chart = samples.plot_bargraph(group_by="tax_id", return_chart=True)
    assert "tax_id" in chart.data.columns
    assert "_tax_id" in chart.data.columns


def test_plot_bargraph_when_passing_tax_name_as_group_by(samples):
    samples._collate_results()
    for sample in samples:
        sample.metadata.custom["tax_name"] = "TAX_NAME"
    chart = samples.plot_bargraph(group_by="tax_name", return_chart=True)
    assert "tax_name" in chart.data.columns
    assert "_tax_name" in chart.data.columns
    assert chart.encoding.tooltip[1].shorthand == "_tax_name"


@pytest.mark.parametrize(
    "method,kwargs",
    [
        ("plot_metadata", {"haxis": "does_not_exist"}),
        ("plot_pca", {"color": "does_not_exist"}),
        ("plot_pca", {"size": "does_not_exist"}),
        ("plot_pca", {"tooltip": "does_not_exist"}),
        ("plot_heatmap", {"tooltip": "does_not_exist"}),
        ("plot_distance", {"tooltip": "does_not_exist"}),
        ("plot_mds", {"tooltip": "does_not_exist"}),
        ("plot_bargraph", {"tooltip": "does_not_exist"}),
    ],
)
def test_plotting_missing_fields(samples, method, kwargs):
    chart = getattr(samples, method)(**kwargs, title="my plot", return_chart=True)
    assert chart.title == "my plot"


def test_interleave_palette_empty_domain():
    assert interleave_palette([]) == []


def test_get_classification_url():
    assert get_classification_url("abc123").endswith("abc123")


@pytest.mark.parametrize(
    "tax_id,expected",
    [
        (None, ""),
        ("", ""),
        ("foo", ""),
        ("2000000000", ""),
        ("123", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=123"),
    ],
)
def test_get_ncbi_taxonomy_browser_url(tax_id, expected):
    assert get_ncbi_taxonomy_browser_url(tax_id) == expected


def test_plot_functional_heatmap(ocx, api_data):
    sample_ids = [
        "543c9c046e3e4e09",
        "66c1531cb0b244f6",
        "37e5151e7bcb4f87",
        "7428cca4a3a04a8e",  # does not have functional results
    ]
    samples = SampleCollection([ocx.Samples.get(x) for x in sample_ids], skip_missing=True)

    assert len(samples) == len(sample_ids)

    with pytest.warns(UserWarning, match="Functional profile not found.*7428cca4a3a04a8e"):
        chart = samples.plot_functional_heatmap(return_chart=True, top_n=3)

    assert len(chart.data.index) == 9  # 3 samples * 3 top_n == 9
    assert set(chart.data["Label"]) == set(x.filename for x in samples[:-1])

    # Defaults to GO metric
    assert all(x.startswith("GO:") for x in chart.data["function_id"])

    assert chart.encoding.x.shorthand == "Label:N"
    assert chart.encoding.y.shorthand == "function_name:N"
    assert {t.shorthand for t in chart.encoding.tooltip} == {
        "Label:N",
        "function_id:N",
        "function_name:N",
        "value:Q",
    }

    # Trying with a different metric
    chart = samples.plot_functional_heatmap(return_chart=True, top_n=3, annotation="eggnog")
    assert all(x.startswith("COG") for x in chart.data["function_id"])


def test_plot_functional_heatmap_only_max_values(ocx, api_data):
    samples = SampleCollection([ocx.Samples.get("543c9c046e3e4e09")])
    chart1 = samples.plot_functional_heatmap(return_chart=True, top_n=2)
    chart2 = samples.plot_functional_heatmap(return_chart=True, top_n=10_000)

    values1 = sorted(list(chart1.data["value"]))
    values2 = sorted(list(chart2.data["value"]))

    assert values1 == values2[len(values2) - 2 :]


def test_plot_functional_heatmap_when_metadata_contains_function_id(ocx, api_data):
    sample_ids = ["543c9c046e3e4e09", "66c1531cb0b244f6", "37e5151e7bcb4f87"]
    samples = SampleCollection([ocx.Samples.get(x) for x in sample_ids])
    for sample in samples:
        sample.metadata.custom["function_id"] = "FUNCTION_ID"

    chart = samples.plot_functional_heatmap(return_chart=True, top_n=3, label="function_id")
    assert set(chart.data["function_name"]) == {
        "[MF] single-stranded DNA endodeoxyribonuclease activity",
        "[CC] phosphopyruvate hydratase complex",
        "[BP] ribosomal large subunit assembly",
    }


@pytest.mark.parametrize(
    ("metadata_key", "escaped_key"),
    [
        ("foo.bar", r"foo\.bar"),
        ("foo[bar]", r"foo\[bar\]"),
    ],
)
def test_plot_bargraph_escape_js_like_fields(samples, metadata_key, escaped_key):
    samples._collate_results(metric="readcount_w_children")
    for i, sample in enumerate(samples):
        sample.metadata.custom[metadata_key] = f"value_{i}"
    chart = samples.plot_bargraph(
        return_chart=True,
        group_by=metadata_key,
        rank="genus",
    ).to_dict()
    assert chart["encoding"]["x"]["field"] == escaped_key


@pytest.mark.parametrize(
    ("metric", "expected_include_taxa_missing_rank"),
    [
        (Metric.Readcount, False),
        (Metric.ReadcountWChildren, True),
        (Metric.Abundance, False),
        (Metric.AbundanceWChildren, True),
    ],
)
def test_plot_bargraph_include_taxa_missing_rank(
    samples, metric, expected_include_taxa_missing_rank
):
    samples._collate_results(metric=metric)
    assert samples._metric == metric

    with (
        mock.patch.object(samples, "to_df") as mock_to_df,
    ):
        samples.plot_bargraph(return_chart=True)

    mock_to_df.assert_called_with(
        rank=Rank.Auto,
        top_n=None,
        threshold=None,
        normalize=None,
        include_taxa_missing_rank=expected_include_taxa_missing_rank,
    )
