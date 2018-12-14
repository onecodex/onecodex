import altair as alt

from onecodex.viz._heatmap import plot_heatmap
from onecodex.viz._pca import plot_pca
from onecodex.viz._distance import plot_distance
from onecodex.viz._metadata import plot_metadata

__all__ = ['plot_heatmap', 'plot_pca', 'plot_distance', 'plot_metadata']


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
alt.renderers.enable('notebook')
