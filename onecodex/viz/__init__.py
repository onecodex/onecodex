import imp
import sys

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
        altair.themes.register("onecodex", onecodex_theme)
        altair.themes.enable("onecodex")

        # Render using `altair_saver` if installed (report environment only, requires node deps)
        if "altair_saver" in altair.renderers.names():
            altair.renderers.enable(
                "altair_saver",
                fmts=["html", "svg"],
                embed_options=VEGAEMBED_OPTIONS,
                vega_cli_options=["--loglevel", "error"],
            )
        else:
            altair.renderers.enable("html", embed_options=VEGAEMBED_OPTIONS)


sys.meta_path = [_AltairImportHook()] + sys.meta_path

__all__ = [
    "VizPCAMixin",
    "VizHeatmapMixin",
    "VizMetadataMixin",
    "VizDistanceMixin",
    "dendrogram",
    "VizBargraphMixin",
]
