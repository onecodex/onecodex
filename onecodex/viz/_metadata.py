import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import matplotlib.dates as mdates

from onecodex.exceptions import OneCodexException
from onecodex.helpers import normalize_classifications, collate_classification_results
from onecodex.distance import alpha_diversity


def plot_metadata(analyses, metadata='created_at', statistic=None, tax_id=None,
                  title=None, label=None, xlabel=None, ylabel=None,
                  field='readcount_w_children', rank='species', normalize=False, **kwargs):
    """Plot by arbitary metadata.

    Note that `rank` only applies if you're calculating a `statistic`.
    Additional **kwargs are passed to Seaborn or Matplotlib plot functions as appropriate.
    """
    if not tax_id and not statistic:
        raise OneCodexException('Please pass a `tax_id` or a `statistic`.')
    elif tax_id and statistic:
        raise OneCodexException('Please pass only a `tax_id` or a `statistic`.')

    sns.set(style=kwargs.pop('style', 'darkgrid'))

    normed_classifications, md = normalize_classifications(analyses, label=label)

    if metadata not in md:
        raise OneCodexException('Selected metadata field `{}` not available'.format(metadata))

    stat = []
    if statistic:
        for analysis in normed_classifications:
            try:
                stat.append(alpha_diversity(analysis, statistic, field=field, rank=rank))
            except ValueError:
                raise NotImplementedError('{} statistic are currently supported.'.format(statistic))
    elif tax_id:
        df, tax_info = collate_classification_results(normed_classifications, field=field, rank=None,
                                                      normalize=normalize)
        stat = df.loc[:, str(tax_id)].values
        if stat.shape[0] == 0:
            raise OneCodexException('No values found for `tax_id` {} and `field` {}.'.format(tax_id, field))

    md['_data'] = stat

    # Try to plot any column with `date` in it as a datetime, or if it's already a datetime
    # Plot numeric types as lmplots
    # Plot categorical, boolean, and objects as boxplot
    if 'date' in metadata.split('_') or pd.api.types.is_datetime64_any_dtype(md[metadata]):
        if not pd.api.types.is_datetime64_any_dtype(md[metadata]):
            md.loc[:, metadata] = md.loc[:, metadata].apply(pd.to_datetime, utc=True)
        fig, ax = plt.subplots()
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
        ax.plot_date(md[metadata].values, md['_data'].values, **kwargs)
        fig.autofmt_xdate()
    elif pd.api.types.is_numeric_dtype(md[metadata]):
        sns.lmplot(x=metadata, y='_data', truncate=True, data=md, **kwargs)
    elif pd.api.types.is_bool_dtype(md[metadata]) or \
            pd.api.types.is_categorical_dtype(md[metadata]) or \
            pd.api.types.is_object_dtype(md[metadata]):
        na = {field: 'N/A' for field in md.columns}
        md.fillna(value=na, inplace=True)
        sns.boxplot(x=metadata, y='_data', data=md, palette='pastel', **kwargs)
    else:
        raise OneCodexException('Unplottable column type for metadata {}'.format(metadata))

    # Labels
    if xlabel is not None:
        plt.gca().set_xlabel(xlabel)
    if ylabel is None:
        ylabel = statistic if statistic else field
    plt.gca().set_ylabel(ylabel)

    _, labels = plt.xticks()
    plt.setp(labels, rotation=90)

    if title:
        plt.suptitle(title)

    plt.show()
