import pandas as pd
import altair as alt

from onecodex.exceptions import OneCodexException


class VizMetadataMixin():
    @staticmethod
    def _box_plot(df, category, quantity, category_type='N',
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

        middle_plot = alt.Chart(df).mark_bar(size=35).encode(
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
            size=35
        ).encode(
            y='median({}):Q'.format(quantity),
            x=alt.X(x_format, axis=alt.Axis(title=xlabel)),
            tooltip='median({}):Q'.format(quantity)
        )

        chart = (lower_plot + middle_plot + upper_plot + middle_tick)

        if title:
            chart = chart.properties(title=title)

        chart.interactive().display()

    def plot_metadata(self, rank='auto',
                      haxis='Label', vaxis='simpson',
                      title=None, xlabel=None, ylabel=None, plot_type='auto'):
        """Plot an arbitrary metadata field versus an arbitrary quantity as a boxplot or scatter plot.

        Specify one of the following for 'haxis' (horizontal axis):
            metadata (string) -- a categorical or numerical metadata field
            tax_id (string) -- taxonomic identifier to be pulled from results
            tax_name (string) -- taxon name to be searched for in results

        In addition to those above, these options are available for 'vaxis' (vertical axis):
            alpha_diversity ('simpson' | 'chao1')
                - 'simpson': calculate alpha diversity using Simpson's Index
                - 'chao1': calculate alpha diversity using Chao1 estimator

        Options for tabulation of classification results:
            rank ('auto' | kingdom' | 'phylum' | 'class' | 'order' | 'family' | 'genus' | 'species')
                - 'auto': choose automatically based on fields
                - 'kingdom' or others: restrict analysis to taxa at this rank

        Options for plotting:
            title (string) -- main title of the plot
            xlabel, ylabel (string) -- axes labels
            plot_type ('auto' | 'boxplot' | 'scatter') -- force a plot type or choose based on data
        """

        if rank is None:
            raise OneCodexException('Please specify a rank or \'auto\' to choose automatically')
        else:
            rank = self._get_auto_rank(rank)

        if plot_type not in ('auto', 'boxplot', 'scatter'):
            raise OneCodexException('Plot type must be one of: auto, boxplot, scatter')

        # alpha diversity is only allowed on vertical axis--horizontal can be magically mapped
        df, magic_fields = self.magic_metadata_fetch([haxis, 'Label'])

        if vaxis in ('simpson', 'chao1'):
            df.loc[:, vaxis] = self.alpha_diversity(vaxis, rank=rank)
            magic_fields[vaxis] = vaxis
        else:
            # if it's not alpha diversity, vertical axis can also be magically mapped
            vert_df, vert_magic_fields = self.magic_metadata_fetch([vaxis])

            # we require the vertical axis to be numerical otherwise plots get weird
            if pd.api.types.is_bool_dtype(vert_df[vert_magic_fields[vaxis]]) or \
               pd.api.types.is_categorical_dtype(vert_df[vert_magic_fields[vaxis]]) or \
               pd.api.types.is_object_dtype(vert_df[vert_magic_fields[vaxis]]) or \
               not pd.api.types.is_numeric_dtype(vert_df[vert_magic_fields[vaxis]]):  # noqa
                raise OneCodexException('Metadata field on vertical axis must be numerical')

            df = pd.concat([df, vert_df], axis=1).dropna(subset=[vert_magic_fields[vaxis]])
            magic_fields.update(vert_magic_fields)

        # plots can look different depending on what the horizontal axis contains
        if pd.api.types.is_datetime64_any_dtype(df[magic_fields[haxis]]):
            category_type = 'T'

            if plot_type == 'auto':
                plot_type = 'boxplot'
        elif 'date' in magic_fields[haxis].split('_'):
            df.loc[:, magic_fields[haxis]] = \
                df.loc[:, magic_fields[haxis]].apply(pd.to_datetime, utc=True)

            category_type = 'T'

            if plot_type == 'auto':
                plot_type = 'boxplot'
        elif pd.api.types.is_bool_dtype(df[magic_fields[haxis]]) or \
             pd.api.types.is_categorical_dtype(df[magic_fields[haxis]]) or \
             pd.api.types.is_object_dtype(df[magic_fields[haxis]]):  # noqa
            df = df.fillna({field: 'N/A' for field in df.columns})

            category_type = 'N'

            if plot_type == 'auto':
                # if data is categorical but there is only one value per sample, scatter plot instead
                if len(df[magic_fields[haxis]].unique()) == len(df[magic_fields[haxis]]):
                    plot_type = 'scatter'
                else:
                    plot_type = 'boxplot'
        elif pd.api.types.is_numeric_dtype(df[magic_fields[haxis]]):
            df = df.dropna(subset=[vaxis])

            category_type = 'O'

            if plot_type == 'auto':
                plot_type = 'scatter'
        else:
            raise OneCodexException('Unplottable column type for horizontal axis ({})'.format(haxis))

        if xlabel is None:
            xlabel = magic_fields[haxis]

        if ylabel is None:
            ylabel = magic_fields[vaxis]

        if plot_type == 'scatter':
            alt.renderers.enable('notebook')

            df = df.reset_index()

            alt_kwargs = dict(
                x=alt.X(magic_fields[haxis], axis=alt.Axis(title=xlabel)),
                y=alt.Y(magic_fields[vaxis], axis=alt.Axis(title=ylabel)),
                tooltip=['Label', '{}:Q'.format(vaxis)],
                href='url:N',
                url='https://app.onecodex.com/classification/' + alt.datum.classification_id
            )

            chart = alt.Chart(df) \
                       .transform_calculate(url=alt_kwargs.pop('url')) \
                       .mark_circle() \
                       .encode(**alt_kwargs)

            if title:
                chart = chart.properties(title=title)

            chart.interactive().display()
        elif plot_type == 'boxplot':
            self._box_plot(
                df,
                magic_fields[haxis],
                magic_fields[vaxis],
                category_type=category_type,
                title=title,
                xlabel=xlabel,
                ylabel=ylabel
            )
