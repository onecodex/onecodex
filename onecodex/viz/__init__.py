import altair as alt
import os

from onecodex.viz._heatmap import VizHeatmapMixin
from onecodex.viz._pca import VizPCAMixin
from onecodex.viz._metadata import VizMetadataMixin
from onecodex.viz._distance import VizDistanceMixin

__all__ = ['VizPCAMixin', 'VizHeatmapMixin', 'VizMetadataMixin', 'VizDistanceMixin']


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


alt.themes.register('onecodex', onecodex_theme)
alt.themes.enable('onecodex')


def onecodex_renderer(spec, **metadata):
    metadata['scale_factor'] = 1.5
    return alt.vegalite.v2.display.png_renderer(spec, **metadata)


alt.renderers.register('onecodex', onecodex_renderer)

if os.getenv('OCX_NBCONVERT', False):
    # if run from notebook service nbconvert, render as a high-res PNG
    alt.renderers.enable('onecodex')
else:
    try:
        # if inside ipython, enable the notebook renderer
        get_ipython()
        alt.renderers.enable('notebook')
    except NameError:
        # otherwise, use altair's default renderer
        alt.renderers.enable('default')
