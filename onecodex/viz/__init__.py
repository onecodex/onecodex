import altair as alt

from onecodex.viz._heatmap import VizHeatmapMixin
from onecodex.viz._pca import VizPCAMixin
from onecodex.viz._primitives import dendrogram
from onecodex.viz._metadata import VizMetadataMixin
from onecodex.viz._distance import VizDistanceMixin
from onecodex.viz._bargraph import VizBargraphMixin


OCX_DARK_GREEN = "#128887"

DEFAULT_PALETTES = {
    "ocx": [
        "#16347B",
        "#0072C7",
        "#01ACEC",
        "#97E9FC",
        "#0A605E",
        "#1DA893",
        "#3DD8BE",
        "#ABEFE2",
        "#37257D",
        "#9C78E0",
        "#CBC0F9",
        "#E3DDFF",
        "#BC5B00",
        "#EB984A",
        "#FCE34D",
        "#FEF2A3",
        "#950303",
        "#DD3A3A",
        "#FF8D8B",
        "#FFD5CB",
        "#771354",
        "#C13A8B",
        "#F28BBF",
        "#F9D9E7",
    ],
    "tableau10": [
        "#4e79a7",
        "#f28e2b",
        "#e15759",
        "#76b7b2",
        "#59a14f",
        "#edc948",
        "#b07aa1",
        "#ff9da7",
        "#9c755f",
        "#bab0ac",
    ],
}

VEGAEMBED_OPTIONS = {
    "mode": "vega-lite",
    "loader": {"target": "_blank", "http": {"credentials": "same-origin"}},
    "logLevel": "error",
}


def onecodex_theme():
    onecodex_palette = ["#ffffcc", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#264153"]

    font_family = "Fira Sans, Helvetica"

    return {
        "config": {
            "range": {"heatmap": list(reversed(onecodex_palette))},
            "area": {"fill": OCX_DARK_GREEN},
            "bar": {"fill": OCX_DARK_GREEN},
            "mark": {"color": OCX_DARK_GREEN},
            "range": {
                "category": DEFAULT_PALETTES["ocx"],
                "heatmap": ["#0A605E", "#1DA893", "#3DD8BE", "#ABEFE2",],
                "ramp": ["#0A605E", "#1DA893", "#3DD8BE", "#ABEFE2",],
            },
            "axis": {
                "labelFont": font_family,
                "labelFontSize": 12,
                "titleFont": font_family,
                "titleFontSize": 12,
                "grid": False,
            },
            "legend": {
                "labelFont": font_family,
                "labelFontSize": 12,
                "titleFont": font_family,
                "titleFontSize": 12,
            },
            "title": {"font": font_family},
            "view": {"width": 400, "height": 400, "strokeWidth": 0},
            "background": "white",
        }
    }


alt.themes.register("onecodex", onecodex_theme)
alt.themes.enable("onecodex")


# Render using `altair_saver` if installed (report environment only, requires node deps)
try:
    import altair_saver  # noqa

    alt.renderers.enable(
        "altair_saver",
        fmts=["html", "svg"],
        embed_options=VEGAEMBED_OPTIONS,
        vega_cli_options=["--loglevel", "error"],
    )
except ImportError:
    alt.renderers.enable("html", embed_options=VEGAEMBED_OPTIONS)

__all__ = [
    "VizPCAMixin",
    "VizHeatmapMixin",
    "VizMetadataMixin",
    "VizDistanceMixin",
    "dendrogram",
    "VizBargraphMixin",
]
