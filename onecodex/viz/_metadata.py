import pandas as pd
import altair as alt
import warnings

from onecodex.exceptions import OneCodexException
from onecodex.helpers import normalize_classifications, collate_classification_results
from onecodex.distance import alpha_diversity


def box_plot(df, category, quantity, category_type='N',
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

    lower_plot = alt.Chart(df).mark_rule().encode(
        y=alt.Y(lower_whisker, axis=alt.Axis(title=ylabel)),
        y2=lower_box,
        x=x_format
    )

    middle_plot = alt.Chart(df).mark_bar(size=35.0).encode(
        y=lower_box,
        y2=upper_box,
        x=x_format
    )

    upper_plot = alt.Chart(df).mark_rule().encode(
        y=upper_whisker,
        y2=upper_box,
        x=x_format
    )

    middle_tick = alt.Chart(df).mark_tick(
        color='black',
        size=35.0
    ).encode(
        y='median({}):Q'.format(quantity),
        x=alt.X(x_format, axis=alt.Axis(title=xlabel), scale=alt.Scale(rangeStep=45)),
        tooltip='median({}):Q'.format(quantity)
    )

    chart = (lower_plot + middle_plot + upper_plot + middle_tick)

    if title:
        chart = chart.properties(title=title)

    chart.interactive().display()


def plot_metadata(analyses, metadata='Label', statistic=None, tax_id=None,
                  title=None, label=None, xlabel=None, ylabel=None,
                  field='readcount_w_children', rank=None, normalize=True,
                  metadata2=None, tax_name=None, boxplot=None, scatter=None):
    """Plot an arbitrary metadata field versus an arbitrary quantity as a boxplot or scatter plot.

    analyses (list) -- list of Samples, Classifications, or Analyses objects to be plotted
    metadata (string) -- metadata field to be plotted on the horizontal axis

    Specify only one of the following options for the vertical axis:
        statistic ('simpson' | 'chao1')
            - 'simpson': calculate alpha diversity using Simpson's Index
            - 'chao1': calculate alpha diversity using Chao1 estimator
        metadata2 (string) -- metadata field to be plotted on the vertical axis
        tax_id (string) -- taxonomic identifier to be pulled from results
        tax_name (string) -- taxon name to be searched for in results

    Options for tabulation of classification results:
        field ('readcount_w_children' | 'readcount' | 'abundance')
            - 'readcount_w_children': total reads of this taxon and all its descendants
            - 'readcount': total reads of this taxon
            - 'abundance': genome size-normalized relative abundances, from shotgun sequencing
        rank ('kingdom' | 'phylum' | 'class' | 'order' | 'family' | 'genus' | 'species')
            - None: include all ranks
            - 'kingdom' or others: restrict analysis to taxa at this rank
        normalize (bool): convert from read counts to relative abundances (each sample sums to 1.0)

    Options for plotting:
        label (string) -- metadata field to label samples with
        title (string) -- main title of the plot
        xlabel, ylabel (string) -- axes labels
        boxplot (bool) -- force output a box plot
        scatter (bool) -- force output a scatter plot
    """

    if sum(map(bool, [statistic, metadata2, tax_id, tax_name])) != 1:
        raise OneCodexException(
            'Please specify exactly one of: statistic, metadata2, tax_id, tax_name'
        )

    if scatter is True and boxplot is True:
        raise OneCodexException('If forcing a plot type, choose only one of boxplot or scatter.')

    normed_classifications, metadata_df = normalize_classifications(analyses, label=label)
    df, tax_info = collate_classification_results(normed_classifications, field=field,
                                                  rank=rank, normalize=normalize)

    metadata_df.index = df.index
    metadata_df['Label'] = metadata_df['_display_name']

    # for displaying what metadata fields are available in the data
    help_fields = ', '.join([str(f) for f in list(metadata_df) + ['classification_id']])

    # figure out what quantity (vertical axis) we're plotting and calculate it
    if statistic in ('simpson', 'chao1'):
        if rank is None:
            rank = 'species'
            warnings.warn('Rank not specified for alpha diversity calculation, using species.')

        metadata_df[statistic] = 0

        for classification in normed_classifications:
            metadata_df.loc[classification.id, statistic] = alpha_diversity(
                classification, statistic, field=field, rank=rank
            )

        vert_axis_field = statistic
    elif statistic is not None:
        raise OneCodexException('Alpha diversity statistic must be one of: simpson, chao1')

    if metadata2:
        if metadata2 not in metadata_df:
            raise OneCodexException(
                'Vertical axis metadata field ({}) not in metadata. Please choose one of the '
                'following fields:\n\n{}'.format(metadata2, help_fields)
            )

        if not pd.api.types.is_numeric_dtype(metadata_df[metadata2]):
            raise OneCodexException(
                'Vertical axis metadata field ({}) must be numerical'.format(metadata2)
            )

        vert_axis_field = metadata2

    if tax_id:
        tax_id = str(tax_id)

        if tax_id not in df:
            if rank:
                raise OneCodexException(
                    'Tax ID {} not found in analyses. Is it in rank={}?'.format(tax_id, rank)
                )
            else:
                raise OneCodexException(
                    'Tax ID {} not found in analyses.'.format(tax_id)
                )

        # if no rank was given, find the rank of this tax_id and re-collate analyses to match
        if rank is None:
            df, tax_info = collate_classification_results(
                normed_classifications,
                field=field,
                rank=tax_info[tax_id]['rank'],
                normalize=normalize
            )

        vert_axis_field = '{} ({})'.format(tax_info[tax_id]['name'], tax_id)
        metadata_df[vert_axis_field] = df[tax_id]

    if tax_name:
        # take the lowest tax_id that maches the search query
        hits = []

        for t in tax_info:
            if tax_name.lower() in tax_info[t]['name'].lower():
                hits.append(t)

        hits = sorted(hits, key=int)

        if not hits:
            if rank:
                raise OneCodexException(
                    'Taxon {} not found in analyses. Is it in rank={}?'.format(tax_name, rank)
                )
            else:
                raise OneCodexException(
                    'Taxon {} not found in analyses.'.format(tax_name)
                )

        # if no rank was given, find the rank of this tax_id and re-collate analyses to match
        if rank is None:
            df, tax_info = collate_classification_results(
                normed_classifications,
                field=field,
                rank=tax_info[hits[0]]['rank'],
                normalize=normalize
            )

        vert_axis_field = '{} ({})'.format(tax_info[hits[0]]['name'], hits[0])
        metadata_df[vert_axis_field] = df[hits[0]]

    # support combination of multiple categorical metadata to plot on x axis
    if isinstance(metadata, list):
        for field in metadata:
            if field not in metadata_df:
                raise OneCodexException(
                    'Horizontal axis metadata field ({}) not in metadata. Please choose one of the '
                    'following fields:\n\n{}'.format(field, help_fields)
                )

            if not (pd.api.types.is_bool_dtype(metadata_df[field]) or  # noqa
                    pd.api.types.is_categorical_dtype(metadata_df[field]) or  # noqa
                    pd.api.types.is_object_dtype(metadata_df[field])):
                raise OneCodexException(
                    'When specifying multiple metadata fields for x, all fields must be categorical'
                )

        composite_field = '_'.join(metadata)

        metadata_df[composite_field] = ''
        metadata_df[composite_field] = metadata_df[composite_field].str.cat(
            [metadata_df[k].astype(str) for k in metadata],
            sep='_'
        ).str.lstrip('_')

        metadata = composite_field
    else:
        if metadata not in metadata_df and metadata != 'classification_id':
            raise OneCodexException(
                'Horizontal axis metadata field ({}) not in metadata. Please choose one of the '
                'following fields:\n\n{}'.format(metadata, help_fields)
            )

    plot_data = pd.DataFrame({
        'classification_id': metadata_df.index,
        metadata: metadata_df[metadata],
        vert_axis_field: metadata_df[vert_axis_field]
    })

    # plots can look different depending on what category_type the data is
    if pd.api.types.is_datetime64_any_dtype(plot_data[metadata]):
        category_type = 'T'

        if not (boxplot or scatter):
            boxplot = True
    elif 'date' in metadata.split('_'):
        plot_data.loc[:, metadata] = plot_data.loc[:, metadata].apply(pd.to_datetime, utc=True)

        category_type = 'T'

        if not (boxplot or scatter):
            boxplot = True
    elif pd.api.types.is_bool_dtype(plot_data[metadata]) or \
         pd.api.types.is_categorical_dtype(plot_data[metadata]) or \
         pd.api.types.is_object_dtype(plot_data[metadata]):  # noqa
        plot_data = plot_data.fillna({field: 'N/A' for field in plot_data.columns})

        category_type = 'N'

        if not (boxplot or scatter):
            # if data is categorical but there is only one value per sample, scatter plot instead
            if len(plot_data[metadata].unique()) == len(plot_data[metadata]):
                scatter = True
            else:
                boxplot = True
    elif pd.api.types.is_numeric_dtype(plot_data[metadata]):
        plot_data = plot_data.dropna(subset=[vert_axis_field])

        category_type = 'O'

        if not (boxplot or scatter):
            scatter = True
    else:
        raise OneCodexException('Unplottable column type for horizontal axis ({})'.format(metadata))

    if xlabel is None:
        xlabel = metadata

    if ylabel is None:
        ylabel = vert_axis_field

    # we can't really tell the difference between data in column 'x' that should be displayed as
    # ordinal or quantitative. should the default behavior for numerical data be scatter?
    if scatter:
        plot_data['Label'] = metadata_df['Label']

        alt_kwargs = dict(
            x=alt.X(metadata, axis=alt.Axis(title=xlabel)),
            y=alt.Y(vert_axis_field, axis=alt.Axis(title=ylabel)),
            tooltip=['Label', '{}:Q'.format(vert_axis_field)],
            href='url:N',
            url='https://app.onecodex.com/classification/' + alt.datum.classification_id
        )

        chart = alt.Chart(plot_data) \
                   .transform_calculate(url=alt_kwargs.pop('url')) \
                   .mark_circle() \
                   .encode(**alt_kwargs)

        if title:
            chart = chart.properties(title=title)

        chart.interactive().display()

    if boxplot:
        box_plot(
            plot_data,
            metadata,
            vert_axis_field,
            category_type=category_type,
            title=title,
            xlabel=xlabel,
            ylabel=ylabel
        )
