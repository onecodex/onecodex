import pytest

pytest.importorskip("pandas")  # noqa

from onecodex.exceptions import OneCodexException


def test_plot_metadata(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    chart = samples.plot_metadata(
        vaxis="simpson", title="my title", xlabel="my xlabel", ylabel="my ylabel", return_chart=True
    )
    assert chart.data["simpson"].round(6).tolist() == [0.788705, 0.845901, 0.372284]
    assert chart.mark == "circle"
    assert chart.title == "my title"
    assert chart.encoding.x.axis.title == "my xlabel"
    assert chart.encoding.y.axis.title == "my ylabel"

    # try time, boolean, and numerical types for x-axis
    chart = samples.plot_metadata(haxis="date_sequenced", vaxis="chao1", return_chart=True)
    assert chart.encoding.x.shorthand == "date_sequenced"

    chart = samples.plot_metadata(haxis="starred", vaxis="chao1", return_chart=True)
    assert chart.encoding.x.shorthand == "starred"

    chart = samples.plot_metadata(haxis="totalige", vaxis="chao1", return_chart=True)
    assert chart.mark == "circle"
    assert chart.encoding.x.shorthand == "totalige"

    # taxid and taxon on vertical axis
    chart = samples.plot_metadata(vaxis=1279, plot_type="boxplot", return_chart=True)
    assert chart.mark.type == "boxplot"

    chart = samples.plot_metadata(vaxis="Staphylococcus", plot_type="scatter", return_chart=True)
    assert chart.mark == "circle"


def test_plot_metadata_exceptions(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    # expect error if rank is None, since that could lead to weird results
    with pytest.raises(OneCodexException) as e:
        samples.plot_metadata(vaxis="simpson", rank=None)
    assert "specify a rank" in str(e.value)

    # vert axis does not exist
    with pytest.raises(OneCodexException):
        samples.plot_metadata(vaxis="does_not_exist")

    # horiz axis does not exist
    with pytest.raises(OneCodexException):
        samples.plot_metadata(haxis="does_not_exist")


def test_plot_pca(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    chart = samples.plot_pca(
        title="my title",
        xlabel="my xlabel",
        ylabel="my ylabel",
        color="geo_loc_name",
        size="totalige",
        org_vectors=3,
        return_chart=True,
        tooltip=["totalige", "vegetables", "Prevotella", 816],
    )

    # one for main plot, one for eigenvectors
    assert len(chart.layer) == 2

    mainplot = chart.layer[0]
    assert mainplot.data["Bacteroides (816)"].round(6).tolist() == [0.356439, 0.213307, 0.787949]
    assert mainplot.data["Prevotella (838)"].round(6).tolist() == [0.000221, 0.001872, 6.7e-05]
    assert mainplot.data["totalige"].tolist() == [62.9, 91.5, 112.0]
    assert mainplot.title == "my title"
    assert mainplot.encoding.color.shorthand == "geo_loc_name"
    assert mainplot.encoding.size.shorthand == "totalige"
    assert mainplot.encoding.x.shorthand == "PC1"
    assert mainplot.encoding.x.axis.title == "my xlabel"
    assert mainplot.encoding.y.shorthand == "PC2"
    assert mainplot.encoding.y.axis.title == "my ylabel"
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


def test_plot_pca_exceptions(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    # expect error if rank is None, since that could lead to weird results
    with pytest.raises(OneCodexException) as e:
        samples.plot_pca(rank=None)
    assert "specify a rank" in str(e.value)

    # must have at least two samples
    with pytest.raises(OneCodexException) as e:
        samples[:1].plot_pca()
    assert "requires 2 or more" in str(e.value)

    # samples must have at least two taxa at this rank
    with pytest.raises(OneCodexException) as e:
        samples.to_df(top_n=1).ocx.plot_pca()
    assert "at least 2 for PCA" in str(e.value)

    # color/size/tooltips with invalid metadata fields or taxids
    for k in ("color", "size", "tooltip"):
        kwargs = {k: "does_not_exist"}
        with pytest.raises(OneCodexException) as e:
            samples.plot_pca(**kwargs)
        assert "not found" in str(e.value)

        kwargs = {k: "487527863"}
        with pytest.raises(OneCodexException) as e:
            samples.plot_pca(**kwargs)
        assert "not found" in str(e.value)


def test_plot_heatmap(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    chart = samples.plot_heatmap(
        top_n=10, title="my title", xlabel="my xlabel", ylabel="my ylabel", return_chart=True
    )
    assert chart.mark == "rect"
    assert chart.title == "my title"
    assert chart.encoding.x.axis.title == "my xlabel"
    assert chart.encoding.y.axis.title == "my ylabel"
    assert len(chart.data["tax_id"].unique()) == 10
    assert chart.data["Reads (Normalized)"].sum().round(6) == 2.613775

    chart = samples.plot_heatmap(threshold=0.1, haxis="eggs", return_chart=True)
    assert chart.mark == "rect"
    assert all(chart.data.groupby("tax_id").max()["Reads (Normalized)"] > 0.1)

    chart = samples.plot_heatmap(top_n=10, threshold=0.01, return_chart=True)
    assert len(chart.data["tax_id"].unique()) == 10
    assert all(chart.data.groupby("tax_id").max()["Reads (Normalized)"] > 0.01)


def test_plot_heatmap_exceptions(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    # expect error if rank is None, since that could lead to weird results
    with pytest.raises(OneCodexException) as e:
        samples.plot_heatmap(rank=None)
    assert "specify a rank" in str(e.value)

    # must have at least two samples
    with pytest.raises(OneCodexException) as e:
        samples[:1].plot_heatmap()
    assert "requires 2 or more" in str(e.value)

    # must specify at least threshold or top_n
    with pytest.raises(OneCodexException) as e:
        samples.plot_heatmap(top_n=None, threshold=None)
    assert "specify at least one of" in str(e.value)

    # tooltip with invalid metadata fields or taxids
    with pytest.raises(OneCodexException) as e:
        samples.plot_heatmap(tooltip="does_not_exist")
    assert "not found" in str(e.value)


def test_plot_distance(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    chart = samples.plot_distance(
        metric="unifrac",
        xlabel="my xlabel",
        ylabel="my ylabel",
        title="my title",
        tooltip="vegetables",
        return_chart=True,
    )
    assert len(chart.hconcat) == 2

    mainplot = chart.hconcat[1]
    assert mainplot.mark == "rect"
    assert mainplot.encoding.x.axis.title == "my xlabel"
    assert mainplot.encoding.y.axis.title == "my ylabel"
    assert sorted([x.shorthand for x in mainplot.encoding.tooltip]) == [
        "1) Label",
        "1) vegetables",
        "2) Label",
        "2) vegetables",
        "Distance:Q",
    ]
    assert chart.hconcat[1].data["Distance"].sum().round(6) == 39.463704


def test_plot_distance_exceptions(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    # expect error if rank is None, since that could lead to weird results
    with pytest.raises(OneCodexException) as e:
        samples.plot_distance(rank=None)
    assert "specify a rank" in str(e.value)

    # only some metrics allowed
    with pytest.raises(OneCodexException) as e:
        samples.plot_distance(metric="simpson")
    assert "must be one of" in str(e.value)

    # need more than one analysis
    with pytest.raises(OneCodexException) as e:
        samples[:1].plot_distance(
            metric="jaccard", xlabel="my xlabel", ylabel="my ylabel", title="my title"
        )
    assert "requires 2 or more" in str(e.value)

    # tooltip with invalid metadata fields or taxids
    with pytest.raises(OneCodexException) as e:
        samples.plot_distance(tooltip="does_not_exist")
    assert "not found" in str(e.value)


def test_plot_mds(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    chart = samples.plot_mds(
        method="pcoa",
        metric="unifrac",
        xlabel="my xlabel",
        ylabel="my ylabel",
        title="my title",
        return_chart=True,
    )
    assert chart.mark == "circle"
    assert chart.title == "my title"
    assert chart.encoding.x.shorthand == "PC1"
    assert chart.encoding.x.axis.title == "my xlabel"
    assert chart.encoding.y.shorthand == "PC2"
    assert chart.encoding.y.axis.title == "my ylabel"
    assert (chart.data["PC1"] * chart.data["PC2"]).sum() == 0.0

    chart = samples.plot_mds(method="smacof", metric="unifrac", return_chart=True)
    assert (chart.data["MDS1"] * chart.data["MDS2"]).sum().round(6) == -0.319449


def test_plot_mds_exceptions(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    # expect error if rank is None, since that could lead to weird results
    with pytest.raises(OneCodexException) as e:
        samples.plot_mds(rank=None)
    assert "specify a rank" in str(e.value)

    # only some metrics allowed
    with pytest.raises(OneCodexException) as e:
        samples.plot_mds(metric="simpson")
    assert "must be one of" in str(e.value)

    # need more than one analysis
    with pytest.raises(OneCodexException) as e:
        samples[:1].plot_mds(
            metric="jaccard", xlabel="my xlabel", ylabel="my ylabel", title="my title"
        )
    assert "requires 2 or more" in str(e.value)

    # tooltip with invalid metadata fields or taxids
    with pytest.raises(OneCodexException) as e:
        samples.plot_mds(tooltip="does_not_exist")
    assert "not found" in str(e.value)


def test_plot_bargraph_arguments(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    # expect error if rank is None, since that could lead to weird results
    with pytest.raises(OneCodexException) as e:
        samples.plot_bargraph(rank=None)
    assert "specify a rank" in str(e.value)

    # tooltip with invalid metadata fields or taxids
    with pytest.raises(OneCodexException) as e:
        samples.plot_bargraph(tooltip="does_not_exist")
    assert "not found" in str(e.value)


def test_plot_bargraph_chart_result(ocx, api_data):
    samples = ocx.Samples.where(project="4b53797444f846c4")

    chart = samples.plot_bargraph(
        rank="phylum",
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
    assert chart.encoding.x.axis.title == "Exemplary Samples"
    assert chart.encoding.y.shorthand == "Reads (Normalized)"
    assert chart.encoding.y.axis.title == "Glorious Abundances"
