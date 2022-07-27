import imp
import sys

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
            "range": {
                "heatmap": list(reversed(onecodex_palette)),
                "categroy": DEFAULT_PALETTES["ocx"],
                "ramp": list(reversed(onecodex_palette)),
            },
            "area": {"fill": OCX_DARK_GREEN},
            "bar": {"fill": OCX_DARK_GREEN},
            "mark": {"color": OCX_DARK_GREEN},
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


def configure_onecodex_theme(altair_module=None):
    """Configure One Codex Altair theme."""
    if altair_module is None:
        try:
            import altair

            altair_module = altair
        except ImportError:
            return  # noop

    altair_module.themes.register("onecodex", onecodex_theme)
    altair_module.themes.enable("onecodex")

    # Render using `altair_saver` if installed (report environment only, requires node deps)
    if "altair_saver" in altair_module.renderers.names():
        # Filter out vega-lite warning about boxplots not yet supporting selection (DEV-4237). This
        # can be removed when vega-lite adds selection support to boxplots:
        #
        # - https://github.com/vega/vega-lite/issues/3702
        # - https://github.com/altair-viz/altair/issues/2232
        def stderr_filter(line):
            """Return ``True`` if stderr line should be displayed."""
            return "Selection not supported for boxplot yet" not in line

        altair_module.renderers.enable(
            "altair_saver",
            fmts=["html", "svg"],
            embed_options=VEGAEMBED_OPTIONS,
            vega_cli_options=["--loglevel", "error"],
            stderr_filter=stderr_filter,
        )
    else:
        altair_module.renderers.enable("html", embed_options=VEGAEMBED_OPTIONS)


# Define an import hook to configure Altair's theme and renderer the first time
# it is imported. Directly importing and configuring Altair in this subpackage
# can slow down the API and CLI. An import hook avoids this performance hit by
# configuring Altair during deferred import in visualization code.
#
# Note: this code is currently Python 2/3 compatible by using the `imp`
# package, which is deprecated in Python 3. Consider using `importlib` if this
# subpackage doesn't need to support Python 2.
#
# Based on: https://stackoverflow.com/a/60352956/3776794
class _AltairImportHook(object):
    def find_module(self, fullname, path=None):
        if fullname != "altair":
            return None
        self.module_info = imp.find_module(fullname, path)
        return self

    def load_module(self, fullname):
        """Load Altair module and configure its theme and renderer."""
        previously_loaded = fullname in sys.modules
        altair = imp.load_module(fullname, *self.module_info)

        if not previously_loaded:
            self._configure_altair(altair)
        return altair

    def _configure_altair(self, altair):
        configure_onecodex_theme(altair)


sys.meta_path = [_AltairImportHook()] + sys.meta_path

__all__ = [
    "VizPCAMixin",
    "VizHeatmapMixin",
    "VizMetadataMixin",
    "VizDistanceMixin",
    "dendrogram",
    "VizBargraphMixin",
]
