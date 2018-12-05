import pandas as pd
import altair as alt
import warnings

from onecodex.exceptions import OneCodexException
from onecodex.helpers import normalize_classifications, collate_classification_results
from onecodex.distance import alpha_diversity


def boxplot(df, category, quantity, category_type='N',
            title=None, xlabel=None, ylabel=None):
    # must be one of Nominal, Ordinal, Time per altair
    if category_type not in ('N', 'O', 'T'):
        raise OneCodexException('If specifying category_type, must be N, O, or T')

    # adapted from https://altair-viz.github.io/gallery/boxplot_max_min.html
    lower_box = 'q1({}):Q'.format(quantity)
    lower_whisker = 'min({}):Q'.format(quantity)
    upper_box = 'q3({}):Q'.format(quantity)
    upper_whisker = 'max({}):Q'.format(quantity)

    if category_type == 'T':
        x_format = 'hoursminutes({}):{}'.format(category, category_type)
    else:
        x_format = '{}:{}'.format(category, category_type)

    if xlabel is None:
        xlabel = category

    if ylabel is None:
        ylabel = quantity

    lower_plot = alt.Chart(df).mark_rule().encode(
        y=alt.Y(lower_whisker, axis=alt.Axis(title=ylabel)),
        y2=lower_box,
        x=x_format
    )

    middle_plot = alt.Chart(df).mark_bar(size=15.0).encode(
        y=lower_box,
        y2=upper_box,
        x=x_format
    )

    upper_plot = alt.Chart(df).mark_rule().encode(
        y=upper_whisker,
        y2=upper_box,
        x=x_format
    )

    # if there is only one sample per category, only middle_tick will be visible. in that case,
    # plot something that is more visible instead of a thin line
    if len(df[category].unique()) == len(df[category]):
        middle_tick = alt.Chart(df).mark_circle(
            size=30.0
        ).encode(
            y='{}:Q'.format(quantity),
            x=alt.X(x_format, axis=alt.Axis(title=xlabel)),
            tooltip='{}:Q'.format(quantity)
        )
    else:
        middle_tick = alt.Chart(df).mark_tick(
            color='black',
            size=15.0
        ).encode(
            y='median({}):Q'.format(quantity),
            x=alt.X(x_format, axis=alt.Axis(title=xlabel)),
            tooltip='median({}):Q'.format(quantity)
        )

    alt.renderers.enable('notebook')

    chart = (lower_plot + middle_plot + upper_plot + middle_tick)

    if title:
        chart = chart.properties(title=title)

    chart.interactive().display()


def plot_metadata(analyses, category='Label', quantity='simpson',
                  title=None, xlabel=None, ylabel=None,
                  field='readcount_w_children', rank=None, normalize=False):
    """Plot an arbitrary metadata field versus an arbitrary quantity as a boxplot.

    analyses (list) -- list of Samples, Classifications, or Analyses objects to be plotted
    category (string) -- metadata field to be plotted on the horizontal axis
    quantity (metadata_field | 'taxid_N' | 'simpson' | 'chao1') -- vertical axis
        - metadata_field: a numerical metadata field
        - 'taxid_N': where N is an arbitrary taxid, report its abundance
        - 'simpson': calculate alpha diversity using Simpson's Index
        - 'chao1': calculate alpha diversity using Chao1 estimator

    field ('readcount_w_children' | 'readcount' | 'abundance')
        - 'readcount_w_children': total reads of this taxon and all its descendants
        - 'readcount': total reads of this taxon
        - 'abundance': genome size-normalized relative abundances, from shotgun sequencing
    rank ('kingdom' | 'phylum' | 'class' | 'order' | 'family' | 'genus' | 'species')
        - 'kingdom' or others: restrict analysis to taxa at this rank
    normalize (bool): convert from read counts to relative abundances (each sample sums to 1.0)

    title (string) -- main title of the plot
    xlabel, ylabel (string) -- axes labels
    """

    if normalize is True and rank is None:
        warnings.warn('Results may be meaningless if rank is None and normalize is True.')

    normed_classifications, metadata = normalize_classifications(analyses)
    df, tax_info = collate_classification_results(normed_classifications, field=field,
                                                  rank=rank, normalize=normalize)

    metadata.index = df.index
    metadata['Label'] = metadata['_display_name']

    # figure out what quantity (vertical axis) we're plotting and calculate it
    if quantity in ('simpson', 'chao1'):
        if rank is None:
            rank = 'species'
            warnings.warn('Rank not specified for alpha diversity calculation, using species.')

        metadata[quantity] = 0

        for classification in normed_classifications:
            metadata.loc[classification.id, quantity] = alpha_diversity(
                classification, quantity, field=field, rank=rank
            )
    elif quantity.startswith('taxid_'):
        taxid = quantity[6:]

        if taxid not in df:
            raise OneCodexException('Tax ID {} not found in analyses'.format(taxid))

        quantity = tax_info[taxid]['name']
        metadata[quantity] = df[taxid]
    elif quantity not in metadata:
        raise OneCodexException(
            'Vertical axis metadata field ({}) not in metadata'.format(quantity)
        )
    elif not pd.api.types.is_numeric_dtype(metadata[quantity]):
        raise OneCodexException(
            'Vertical axis metadata field ({}) must be numerical'.format(quantity)
        )

    metadata = metadata.reset_index()

    if category not in metadata:
        raise OneCodexException(
            'Horizontal axis metadata field ({}) not in metadata'.format(category)
        )

    plot_data = pd.DataFrame(data={category: metadata[category], quantity: metadata[quantity]})

    # we're going to use boxplot() no matter what, but we can pass different values
    # for category_type which will determine the plot we get
    if pd.api.types.is_datetime64_any_dtype(plot_data[category]):
        category_type = 'T'
    elif 'date' in category.split('_'):
        plot_data.loc[:, category] = plot_data.loc[:, category].apply(pd.to_datetime, utc=True)
        category_type = 'T'
    elif pd.api.types.is_bool_dtype(plot_data[category]) or \
         pd.api.types.is_categorical_dtype(plot_data[category]) or \
         pd.api.types.is_object_dtype(plot_data[category]):  # noqa
        plot_data = plot_data.fillna({field: 'N/A' for field in plot_data.columns})
        category_type = 'N'
    elif pd.api.types.is_numeric_dtype(plot_data[category]):
        plot_data = plot_data.dropna(subset=[quantity])
        category_type = 'O'
    else:
        raise OneCodexException('Unplottable column type for category {}'.format(category))

    boxplot(
        plot_data,
        category,
        quantity,
        category_type=category_type,
        title=title,
        xlabel=xlabel,
        ylabel=ylabel
    )
