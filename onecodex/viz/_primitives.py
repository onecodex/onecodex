import altair as alt
import pandas as pd

from onecodex.exceptions import OneCodexException


def sort_helper(sort, values):
    sort_order = None
    if callable(sort):
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
