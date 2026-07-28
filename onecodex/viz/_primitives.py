from __future__ import annotations

import os
from itertools import chain
from collections.abc import Iterable
from math import ceil
from typing import TYPE_CHECKING

from onecodex.exceptions import OneCodexException

if TYPE_CHECKING:
    import altair as alt


# We use custom IDs in some databases and we don't want to link them to NCBI. They are always after 2e9
CUSTOM_TAX_ID_START = 2e9


def get_classification_url(classification_id):
    return f"{os.environ.get('ONE_CODEX_API_BASE', 'https://app.onecodex.com')}/classification/{classification_id}"


def get_ncbi_taxonomy_browser_url(tax_id):
    try:
        tax_id = int(tax_id)
    except (TypeError, ValueError):  # "Other", "No genus", etc.
        return ""

    if tax_id < CUSTOM_TAX_ID_START:
        return f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={tax_id}"
    return ""


def sort_helper(sort, values):
    """Return a sorted list of values for the Altair chart axes."""
    sort_order = None

    if callable(sort):
        values = list(set(values))
        sort_order = sort(values)
    elif isinstance(sort, list):
        if set(sort) != set(values):
            raise OneCodexException("sort_x must have the same items as your dataset.")
        sort_order = sort
    elif sort:
        raise OneCodexException(
            "Please pass either a sorted list of values matching the axis labels \
            or a function that returns a sorted list of labels"
        )

    return sort_order


def prepare_props(title=None, height=None, width=None):
    """Prepare key plotting kwargs for passing to Altair, which d/n like None values."""
    props = {}
    if title:
        props["title"] = title  # None gets rendered as `None`
    if height:
        props["height"] = height  # None violates Vega JSON Schema
    if width:
        props["width"] = width  # None violates Vega JSON Schema
    return props


def interleave_palette(domain, palette="ocx"):
    from onecodex.viz import DEFAULT_PALETTES

    if palette in DEFAULT_PALETTES:
        colors = DEFAULT_PALETTES[palette]
    elif isinstance(palette, list):
        colors = palette
    else:
        raise OneCodexException("A valid palette name or list of colors must be passed")

    n_rows = len(set(domain))
    if n_rows == 0:
        return []

    # We do some shuffling to optimize the range of colours with our own palette
    if palette == "ocx":
        hues, shades = 6, 4

        # Calculate how many shades to show of each hue
        period = min(ceil(n_rows / hues), shades)

        # Save the darkest hue for last
        offset = 0
        if period < shades:
            offset = 1

        # Generate a sub-palette with each hue, at each shade
        sub_palettes = [colors[ix::shades] for ix in range(offset, period + offset)]

        # Interleave the sub-palettes so the individuals colors are ordered by hue, then shade
        colors = list(chain.from_iterable(zip(*sub_palettes)))

    # Repeat the palette to extend it to the length of the domain
    extended_palette = colors * (n_rows // len(colors)) + colors[: n_rows % len(colors)]

    return extended_palette


def dendrogram(tree):
    """Plot a simple square dendrogram using Altair.

    Parameters
    ----------
    tree : `dict` returned by `scipy.cluster.hierarchy.dendrogram`
    Contains, at a minimum, 'icoord', 'dcoord', and 'leaves' keys. Scipy does all the work of
    determining where the lines in the tree should go. All we have to do is draw them.

    Returns
    -------
    `altair.Chart`
    """
    # Deferred imports
    import altair as alt
    import pandas as pd

    plot_data = {
        "x": [],
        "y": [],
        "o": [],  # order these points should be connected in
        "b": [],  # one number per branch
    }

    for idx, (i, d) in enumerate(zip(tree["icoord"], tree["dcoord"])):
        plot_data["x"].extend(map(lambda x: -x, d))
        plot_data["y"].extend(map(lambda x: -x, i))
        plot_data["o"].extend([0, 1, 2, 3])
        plot_data["b"].extend([idx] * 4)

    plot_data = pd.DataFrame(plot_data)

    chart = (
        alt.Chart(plot_data, width=100, height=15 * len(tree["leaves"]) - 7.5)
        .mark_line(point=False, opacity=0.5)
        .encode(
            x=alt.X("x", axis=None),
            y=alt.Y("y", axis=None, scale=alt.Scale(zero=True, nice=False)),
            order="o",
            color=alt.Color(
                "b:N",
                scale=alt.Scale(domain=list(range(idx + 1)), range=["black"] * (idx + 1)),
                legend=None,
            ),
        )
    )

    return chart


def get_unique_column(preferred_name, existing_columns):
    """Generate a unique column name.

    Parameters
    ----------
    preferred_name : `str`
        Starting point. Returns this if it does not collide with `existing_columns`.
    existing_columns : `list`

    Returns
    -------
    `str`
    """
    result = preferred_name
    while result in existing_columns:
        result = "_" + result
    return result


# Altair shorthand type codes (trailing ":<code>", e.g. "foo:N")
_ALTAIR_TYPECODES = frozenset(
    ["O", "N", "Q", "T", "G", "ordinal", "nominal", "quantitative", "temporal", "geojson"]
)


def _escape_field(field: str) -> str:
    """Escape nested object/array access path characters for Vega-Lite's field accessor."""
    return field.replace(".", r"\.").replace("[", r"\[").replace("]", r"\]")


def _field_name_has_colon(shorthand: str, columns: set[str]) -> bool:
    """Report whether the shorthand's field name contains a colon Altair would misread as a type."""
    if ":" not in shorthand:
        return False

    if shorthand in columns:
        # column literally named e.g. "foo:N"
        return True

    field, sep, maybe_type = shorthand.rpartition(":")
    real_field = field if (sep and maybe_type in _ALTAIR_TYPECODES) else shorthand
    return ":" in real_field


def escape_chart_fields(chart: alt.Chart) -> None:
    """Escape characters in field names that Altair/Vega-Lite would otherwise treat as syntax.

    WARNING: it modifies `chart` in-place.

    Parameters
    ----------
    chart : `altair.Chart`

    """
    import altair as alt
    import pandas as pd
    from altair.utils import parse_shorthand

    data = chart.data if isinstance(getattr(chart, "data", None), pd.DataFrame) else None
    columns = set(data.columns) if data is not None else set()

    def _to_explicit_field(schema_item: alt.VegaLiteSchema, shorthand: str) -> None:
        schema_item._kwds["shorthand"] = alt.Undefined
        schema_item._kwds["field"] = _escape_field(shorthand)
        if shorthand not in columns:
            # Not in the data, let Vega-Lite infer the type at render
            return

        # Infer type/sort if they aren't explicitly declared
        infer_type = schema_item._kwds.get("type") is alt.Undefined
        infer_sort = schema_item._kwds.get("sort") is alt.Undefined
        if not infer_type and not infer_sort:
            return

        # Infer type/sort using Altair
        attrs = parse_shorthand({"field": shorthand.replace(":", r"\:")}, data)
        if infer_type and "type" in attrs:
            schema_item._kwds["type"] = attrs["type"]
        if infer_sort and "sort" in attrs:
            schema_item._kwds["sort"] = attrs["sort"]

    def _escape_iter(schema_item: alt.VegaLiteSchema) -> None:
        for key, val in schema_item._kwds.items():
            if isinstance(val, alt.VegaLiteSchema):
                _escape_iter(val)
            elif isinstance(val, Iterable) and not isinstance(val, str):
                for v in val:
                    if isinstance(v, alt.VegaLiteSchema):
                        _escape_iter(v)

            elif key == "shorthand" and isinstance(val, str):
                if _field_name_has_colon(val, columns):
                    _to_explicit_field(schema_item, val)
                else:
                    schema_item._kwds[key] = _escape_field(val)

    if isinstance(chart.encoding, alt.VegaLiteSchema):
        _escape_iter(chart.encoding)
