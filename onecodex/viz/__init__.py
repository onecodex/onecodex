import altair as alt
from functools import partial
import json
import uuid

from onecodex.exceptions import OneCodexException
from onecodex.viz._heatmap import VizHeatmapMixin
from onecodex.viz._pca import VizPCAMixin
from onecodex.viz._primitives import dendrogram, boxplot
from onecodex.viz._metadata import VizMetadataMixin
from onecodex.viz._distance import VizDistanceMixin
from onecodex.viz._bargraph import VizBargraphMixin


# # public CDN
# ONE_CODEX_VEGA_CDN = 'https://cdn.jsdelivr.net/npm/'
ONE_CODEX_VEGA_CDN = "https://static.onecodex.com/cdn/"


def onecodex_theme():
    onecodex_palette = ["#ffffcc", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#264153"]

    return {
        "config": {
            "range": {"heatmap": list(reversed(onecodex_palette))},
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


def onecodex_renderer(spec, png=None, svg=None, save_json=None, **metadata):
    if png is None or svg is None or save_json is None:
        raise OneCodexException("One Codex Vega renderer is not properly configured")

    # see https://github.com/vega/vega-embed/issues/8
    vega_js = """
// store the vega spec in an element unique to this figure. can't rely on `element` being the
// same inside require() because it waits until all the modules load, at which point the next
// figure's output area will have been rendered.
var out_target = document.createElement("div");
out_target.id = "{id}";
out_target.vega_spec = {spec};
out_target.className = "vega-embed";

var style = document.createElement("style");
style.textContent = [
".vega-embed .error p {{",
"  color: firebrick;",
"  font-size: 14px;",
"}}",
].join("\\n");

if (element[0]) {{
  // jupyter notebook uses jquery, which returns a list of elements
  element[0].appendChild(out_target);
  element[0].appendChild(style);
}} else {{
  // jupyter lab returns the DOM element directly
  element.appendChild(out_target);
  element.appendChild(style);

  // requirejs is not available in jupyter lab
  var require_script = document.getElementById("one-codex-requirejs");

  if (!require_script) {{
    require_script = document.createElement("script");
    require_script.id = "one-codex-requirejs";
    require_script.src = "{cdn}require.min.js";
    document.head.appendChild(require_script);
  }}
}}

var output_area = this;

requirejs.config({{
  baseUrl: '{cdn}',
  paths: {{
    "vega-embed":  "vega-embed@4?noext",
    "vega-lib": "vega-lib@4?noext",
    "vega-lite": "vega-lite@3?noext",
    "vega": "vega@5?noext"
  }}
}});

function showError(el, error) {{
  el.innerHTML = `<div class="error">
    <p>JavaScript Error: ${{error.message}}</p>
    <p>This usually means there's a typo in your chart specification.
    See the JavaScript console for the full traceback.</p>
  </div>`;

  throw error;
}}

require(["vega-embed"], function(vegaEmbed) {{
  out_target = document.getElementById("{id}");

  vegaEmbed(out_target, out_target.vega_spec, {{
    defaultStyle: true,
    loader: {{ target: "_blank", http: {{ credentials: "same-origin" }} }},
    mode: "vega-lite"
  }}){image_save_js}.catch(error => showError(out_target, error));
}}, function (err) {{
  if (err.requireType !== "scripterror") {{
    throw(err);
  }}
}});
    """

    image_save_js = """
.then(result => {{
    const imageData = result.view.{image_save_func}.then(imageData => {{
      if (output_area !== undefined && output_area.outputs !== undefined) {{
        // figure out which output cell belongs to this render block. there may be
        // multiple jupyter-vega cells per input cell, but only one will match our id
        for (const cell_num in output_area.outputs) {{
          let cell = output_area.outputs[cell_num];
          if (cell.metadata && cell.metadata["jupyter-vega"] === "#{id}") {{
            output_area.outputs[cell_num]["data"]["{image_mime_type}"] = imageData;
          }}
        }}
      }}
    }}).catch(error => showError(out_target, error));
  }})
"""
    element_uuid = uuid.uuid4().hex

    # CSS selector ids cannot start with numbers. if ours does, rename it
    if element_uuid[0] in map(str, range(10)):
        element_uuid = "a" + element_uuid[1:]

    # figure out what to call to save SVG/PNGs
    if svg:
        image_save_func = "toSVG()"
        image_mime_type = "image/svg+xml"

    if png:
        image_save_func = "toImageURL('png')"
        image_mime_type = "image/png"

    if svg or png:
        image_save_js = image_save_js.format(
            id=element_uuid, image_save_func=image_save_func, image_mime_type=image_mime_type
        )
    else:
        image_save_js = ""

    ret_metadata = {"jupyter-vega": "#{0}".format(element_uuid)}
    ret_dict = {
        "application/javascript": vega_js.format(
            id=element_uuid,
            spec=json.dumps(spec),
            cdn=ONE_CODEX_VEGA_CDN,
            image_save_js=image_save_js,
        )
    }

    if save_json:
        ret_dict["application/vnd.vegalite.v2+json"] = json.dumps(spec)

    return (ret_dict, ret_metadata)


def renderer_settings(svg_or_png=None, save_json=True, enable=True):
    """Change behavior of Vega/Altair renderer in IPython notebook for this session.

    Parameters
    ----------
    svg_or_png : `str` in {"png", "svg"} or `False` to disable saving images
        Save rendered image in PNG or SVG format in an output cell. Defaults to "svg"
    save_json : `bool`
        Store Altair-generated JSON in output cell.
    enable : `bool`
        If True, after updating renderer settings, will enable the renderer. If False, user must
        call `altair.renderers.enable("onecodex")` before changes will take effect.
    """
    svg = png = False

    if svg_or_png is None or svg_or_png == "svg":
        svg = True
    elif svg_or_png == "png":
        png = True
    elif svg_or_png is False:
        pass
    else:
        raise OneCodexException("svg_or_png kwarg must be one of: png, svg")

    renderer = partial(onecodex_renderer, svg=svg, png=png, save_json=save_json)
    alt.renderers.register("onecodex", renderer)

    if enable:
        alt.renderers.enable("onecodex")


alt.themes.register("onecodex", onecodex_theme)
alt.themes.enable("onecodex")
renderer_settings(enable=False)


try:
    # if inside ipython, enable the notebook renderer
    get_ipython()
    alt.renderers.enable("onecodex")
except NameError:
    # otherwise, use altair's default renderer
    alt.renderers.enable("default")


__all__ = [
    "VizPCAMixin",
    "VizHeatmapMixin",
    "VizMetadataMixin",
    "VizDistanceMixin",
    "dendrogram",
    "boxplot",
    "VizBargraphMixin",
]
