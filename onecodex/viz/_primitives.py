import altair as alt
import pandas as pd

from onecodex.exceptions import OneCodexException


def boxplot(df, category, quantity, category_type="N", title=None, xlabel=None, ylabel=None):
    """Plot a simple boxplot using Altair.

    Parameters
    ----------
    df : `pandas.DataFrame`
        Contains columns matching 'category' and 'quantity' labels, at a minimum.
    category : `string`
        The name of the column in df used to group values on the horizontal axis.
    quantity : `string`
        The name of the columm in df of values to plot on the vertical axis. Must be numerical.
    category_type : {'N', 'O', 'T'}, optional
        Nominal, ordinal, or time values can be used as categories. Quantitative (Q) values look weird.
    title : `string`, optional
        Text label at the top of the plot.
    xlabel : `string`, optional
        Text label along the horizontal axis.
    ylabel : `string`, optional
        Text label along the vertical axis.

    Returns
    -------
    `altair.Chart`
    """
    # must be one of Nominal, Ordinal, Time per altair
    if category_type not in ("N", "O", "T"):
        raise OneCodexException("If specifying category_type, must be N, O, or T")

    # adapted from https://altair-viz.github.io/gallery/boxplot_max_min.html
    lower_box = "q1({}):Q".format(quantity)
    lower_whisker = "min({}):Q".format(quantity)
    upper_box = "q3({}):Q".format(quantity)
    upper_whisker = "max({}):Q".format(quantity)

    if category_type == "T":
        x_format = "hoursminutes({}):{}".format(category, category_type)
    else:
        x_format = "{}:{}".format(category, category_type)

    lower_plot = (
        alt.Chart(df)
        .mark_rule()
        .encode(y=alt.Y(lower_whisker, axis=alt.Axis(title=ylabel)), y2=lower_box, x=x_format)
    )

    middle_plot = alt.Chart(df).mark_bar(size=35).encode(y=lower_box, y2=upper_box, x=x_format)

    upper_plot = alt.Chart(df).mark_rule().encode(y=upper_whisker, y2=upper_box, x=x_format)

    middle_tick = (
        alt.Chart(df)
        .mark_tick(color="black", size=35)
        .encode(
            y="median({}):Q".format(quantity),
            x=alt.X(x_format, axis=alt.Axis(title=xlabel), scale=alt.Scale(rangeStep=45)),
            tooltip="median({}):Q".format(quantity),
        )
    )

    chart = lower_plot + middle_plot + upper_plot + middle_tick

    if title:
        chart = chart.properties(title=title)

    return chart


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
