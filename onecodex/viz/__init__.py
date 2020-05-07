import altair as alt

from onecodex.viz._heatmap import VizHeatmapMixin
from onecodex.viz._pca import VizPCAMixin
from onecodex.viz._primitives import dendrogram
from onecodex.viz._metadata import VizMetadataMixin
from onecodex.viz._distance import VizDistanceMixin
from onecodex.viz._bargraph import VizBargraphMixin


OCX_DARK_GREEN = "#128887"
# OCX_VEGA_CDN = "https://static.onecodex.com/cdn"


VEGAEMBED_OPTIONS = {
    "mode": "vega-lite",
    "loader": {"target": "_blank", "http": {"credentials": "same-origin"}},
    "logLevel": "error",
}


def onecodex_theme():
    onecodex_palette = ["#ffffcc", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#264153"]

    return {
        "config": {
            "range": {"heatmap": list(reversed(onecodex_palette))},
            "area": {"fill": OCX_DARK_GREEN},
            "bar": {"fill": OCX_DARK_GREEN},
            "axis": {
                "labelFont": "Helvetica",
                "labelFontSize": 12,
                "titleFont": "Helvetica",
                "titleFontSize": 12,
                "grid": False,
            },
            "legend": {
                "labelFont": "Helvetica",
                "labelFontSize": 12,
                "titleFont": "Helvetica",
                "titleFontSize": 12,
            },
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
