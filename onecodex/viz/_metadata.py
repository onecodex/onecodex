import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import matplotlib.dates as mdates

from onecodex.exceptions import OneCodexException
from onecodex.helpers import normalize_classifications, collate_classification_results
from onecodex.distance import alpha_diversity


def plot_metadata(analyses, metadata='created_at', statistic=None, tax_id=None,
                  title=None, label=None, xlabel=None, ylabel=None,
                  field='readcount_w_children', rank='species'):
    # metadata -> string (metadata field in default *or* custom metadata) -> x axis?
    # metadata -> (lambda?) # TODO
    # ONE of:
    # tax_id -> plot abundance on y axis
    # statistic -> plot statistic on y axis
    # statistic (lambda?) # TODO
    if not tax_id and not statistic:
        raise OneCodexException('Please pass a `tax_id` or a `statistic`.')
    elif tax_id and statistic:
        raise OneCodexException('Please pass only a `tax_id` or a `statistic`.')

    sns.set(style="whitegrid")

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
        df = collate_classification_results(normed_classifications, field=field, rank=rank)
        stat = df.loc[:, str(tax_id)].values

    md['_data'] = stat

    # Try to plot any column with `date` in it as a datetime, or if it's already a datetime
    # Plot numeric types as lmplots
    # Plot categorical, boolean, and objects as boxplot
    if 'date' in metadata.split('_') or pd.core.dtypes.common.is_datetime64_any_dtype(md[metadata]):
        if not pd.core.dtypes.common.is_datetime64_any_dtype(md[metadata]):
            md.loc[:, metadata] = md.loc[:, metadata].apply(pd.to_datetime, utc=True)
        fig, ax = plt.subplots()
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
        ax.plot_date(md[metadata].values, md['_data'].values)
        fig.autofmt_xdate()
    elif pd.core.dtypes.common.is_numeric_dtype(md[metadata]):
        sns.lmplot(x=metadata, y='_data', truncate=True, data=md)
    elif pd.core.dtypes.common.is_bool_dtype(md[metadata]) or \
            pd.core.dtypes.common.is_categorical_dtype(md[metadata]) or \
            pd.core.dtypes.common.is_object_dtype(md[metadata]):
        md.fillna(value='N/A', inplace=True)
        sns.boxplot(x=metadata, y='_data', data=md, palette='pastel')
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
