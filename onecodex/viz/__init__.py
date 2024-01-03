import importlib.abc
import importlib.util
import sys

from onecodex.viz._heatmap import VizHeatmapMixin
from onecodex.viz._pca import VizPCAMixin
from onecodex.viz._primitives import dendrogram
from onecodex.viz._metadata import VizMetadataMixin
from onecodex.viz._distance import VizDistanceMixin
from onecodex.viz._bargraph import VizBargraphMixin
from onecodex.viz._functional import VizFunctionalHeatmapMixin


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
                "category": DEFAULT_PALETTES["ocx"],
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
        import functools
        import shutil
        import altair_saver.savers._node
        from altair_saver._utils import check_output_with_stderr

        # Change `npm bin` to `npm root` for compatibility with npm >=9 (also backwards compatible
        # with npm 8). We are monkeypatching because altair-saver appears to be an unmaintained
        # project. We are pinned at v0.5.0 so this hack should be safe.
        #
        # Function copied and modified from:
        # https://github.com/altair-viz/altair_saver/blob/v0.5.0/altair_saver/savers/_node.py#L15-L24
        #
        # Applied the fix from:
        # https://github.com/altair-viz/altair_saver/pull/116
        #
        # altair-saver is BSD-3-Clause:
        # https://github.com/altair-viz/altair_saver/blob/v0.5.0/LICENSE
        #
        # altair-saver license is included with this software in `licenses/altair-saver.txt`
        @functools.lru_cache(2)
        def npm_bin(global_: bool) -> str:
            """Locate the npm binary directory."""
            npm = shutil.which("npm")
            if not npm:
                raise altair_saver.savers._node.ExecutableNotFound("npm")
            cmd = [npm, "root"]
            if global_:
                cmd.append("--global")
            return check_output_with_stderr(cmd).decode().strip()

        altair_saver.savers._node.npm_bin = npm_bin

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
# Some inspiration taken from:
# https://github.com/GrahamDumpleton/wrapt/blob/master/src/wrapt/importer.py
class _AltairFinder(importlib.abc.MetaPathFinder):
    def __init__(self):
        self.in_progress = False

    def find_spec(self, fullname, path, target=None):
        # `importlib.util.find_spec()` searches `sys.meta_path` for a module spec, which includes
        # `self`. Use `in_progress` flag to avoid infinite recursion.
        if fullname != "altair" or self.in_progress:
            return None

        self.in_progress = True
        try:
            spec = importlib.util.find_spec(fullname)
            spec.loader = _AltairLoader(spec.loader)
            return spec
        finally:
            self.in_progress = False


class _AltairLoader(importlib.abc.Loader):
    def __init__(self, loader):
        self.loader = loader

    def create_module(self, spec):
        return None  # use default module creation semantics

    def exec_module(self, module):
        self.loader.exec_module(module)
        configure_onecodex_theme(module)


sys.meta_path = [_AltairFinder()] + sys.meta_path

__all__ = [
    "VizPCAMixin",
    "VizHeatmapMixin",
    "VizMetadataMixin",
    "VizDistanceMixin",
    "dendrogram",
    "VizBargraphMixin",
    "VizFunctionalHeatmapMixin",
]
