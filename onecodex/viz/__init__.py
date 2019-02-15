import altair as alt
import json
import uuid

from onecodex.viz._heatmap import VizHeatmapMixin
from onecodex.viz._pca import VizPCAMixin
from onecodex.viz._primitives import dendrogram, boxplot
from onecodex.viz._metadata import VizMetadataMixin
from onecodex.viz._distance import VizDistanceMixin
from onecodex.viz._bargraph import VizBargraphMixin


# # public CDN
# ONE_CODEX_VEGA_CDN = 'https://cdn.jsdelivr.net/npm/'
ONE_CODEX_VEGA_CDN = 'https://static.onecodex.com/cdn/'


def onecodex_theme():
    onecodex_palette = [
        '#ffffcc',
        '#c7e9b4',
        '#7fcdbb',
        '#41b6c4',
        '#2c7fb8',
        '#264153'
    ]

    return {
        'config': {
            'range': {
                'heatmap': list(reversed(onecodex_palette))
            },
            'axis': {
                'labelFont': 'Helvetica',
                'labelFontSize': 12,
                'titleFont': 'Helvetica',
                'titleFontSize': 12,
                'grid': False
            },
            'legend': {
                'labelFont': 'Helvetica',
                'labelFontSize': 12,
                'titleFont': 'Helvetica',
                'titleFontSize': 12
            },
            'view': {
                'width': 400,
                'height': 400,
                'strokeWidth': 0
            },
            'background': 'white'
        }
    }


def onecodex_renderer(spec, **metadata):
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
    "vega-embed":  "vega-embed@3?noext",
    "vega-lib": "vega-lib?noext",
    "vega-lite": "vega-lite@2?noext",
    "vega": "vega@3?noext"
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
  }})
    .then(result => {{
      const imageData = result.view
        .toSVG()
        .then(imageData => {{
          if (output_area !== undefined && output_area.outputs !== undefined) {{
            // figure out which output cell belongs to this render block. there may be
            // multiple jupyter-vega cells per input cell, but only one will match our id
            for (const cell_num in output_area.outputs) {{
              let cell = output_area.outputs[cell_num];
              if (cell.metadata && cell.metadata["jupyter-vega"] === "#{id}") {{
                output_area.outputs[cell_num]["data"]["image/svg+xml"] = imageData;
              }}
            }}
          }}
        }})
        .catch(error => showError(out_target, error));
    }})
    .catch(error => showError(out_target, error));
}}, function (err) {{
  if (err.requireType !== "scripterror") {{
    throw(err);
  }}
}});
    """
    el_uuid = uuid.uuid4().hex

    # CSS selector ids cannot start with numbers. if ours does, rename it
    if el_uuid[0] in map(str, range(10)):
        el_uuid = 'a' + el_uuid[1:]

    return (
        {
            'application/javascript': vega_js.format(
                id=el_uuid, spec=json.dumps(spec), cdn=ONE_CODEX_VEGA_CDN
            ),
            'application/vnd.vegalite.v2+json': json.dumps(spec)
        },
        {'jupyter-vega': '#{0}'.format(el_uuid)}
    )


alt.themes.register('onecodex', onecodex_theme)
alt.themes.enable('onecodex')
alt.renderers.register('onecodex', onecodex_renderer)


try:
    # if inside ipython, enable the notebook renderer
    get_ipython()
    alt.renderers.enable('onecodex')
except NameError:
    # otherwise, use altair's default renderer
    alt.renderers.enable('default')


__all__ = ['VizPCAMixin', 'VizHeatmapMixin', 'VizMetadataMixin', 'VizDistanceMixin',
           'dendrogram', 'boxplot', 'VizBargraphMixin']
