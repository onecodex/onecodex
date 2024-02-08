import pytest
import mock

pytest.importorskip("pandas")  # noqa

import pandas as pd

from onecodex.exceptions import StatsException, PlottingException


def test_differential_abundance_not_enough_groups(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")
    samples.metadata["wheat"] = "one value"

    with pytest.raises(StatsException, match="`group_by` must have at least two groups"):
        samples.differential_abundance(group_by="wheat")


def test_differential_abundance_no_group_variance(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")
    samples.metadata["wheat"] = ["all", "unique", "values"]

    with pytest.raises(StatsException, match="`group_by` contains only one value"):
        samples.differential_abundance(group_by="wheat")


def test_differential_abundance_one_group_var(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")
    results = samples.differential_abundance(group_by="wheat")

    assert results.ancom_df.index.name == "tax_id"
    assert "W" in results.ancom_df.columns
    assert "Reject null hypothesis" in results.ancom_df.columns
    assert results.ancom_df["Reject null hypothesis"].any()
    assert results.percentile_df.index.name == "tax_id"

    chart = results.plot(return_chart=True)
    assert chart.facet.shorthand == "Significantly Different Taxa"
    assert "wheat" in chart.data.columns


def test_differential_abundance_multiple_group_vars(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")
    samples.metadata["tofu"] = ["soft", "extra firm", "soft"]
    results = samples.differential_abundance(group_by=("wheat", "tofu"))

    assert results.ancom_df.index.name == "tax_id"
    assert "W" in results.ancom_df.columns
    assert "Reject null hypothesis" in results.ancom_df.columns
    assert results.ancom_df["Reject null hypothesis"].any()
    assert results.percentile_df.index.name == "tax_id"

    chart = results.plot(return_chart=True)
    assert chart.facet.shorthand == "Significantly Different Taxa"
    assert "wheat_tofu" in chart.data.columns


def test_differential_abundance_no_significant_taxa(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")
    ancom_df = pd.DataFrame(
        [[0, False], [5, False], [1, False]], columns=["W", "Reject null hypothesis"]
    )
    percentile_df = pd.DataFrame()

    with mock.patch("skbio.stats.composition.ancom", return_value=(ancom_df, percentile_df)):
        results = samples.differential_abundance(group_by="wheat")

    assert not results.ancom_df["Reject null hypothesis"].any()

    with pytest.raises(PlottingException, match="No significantly different taxa"):
        results.plot()


def test_differential_abundance_plot_props(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")
    results = samples.differential_abundance(group_by="wheat")
    chart = results.plot(title="New Title", width=100, height=100, return_chart=True)

    assert chart.title == "New Title"
    assert chart.width == 100
    assert chart.height == 100
