import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib.dates as mdates

from onecodex.helpers import normalize_analyses, collate_analysis_results
from onecodex.distance import alpha_diversity


def plot_metadata(analyses, title=None, metadata='created_at', statistic=None,
                  tax_id=None, field='readcount_w_children', rank='species'):
    # metadata -> string (metadata field in default *or* custom metadata) -> x axis?
    # metadata -> (lambda?) # TODO
    # ONE of:
    # tax_id -> plot abundance on y axis
    # statistic -> plot statistic on y axis
    # statistic (lambda?) # TODO
    assert tax_id or statistic
    assert not (tax_id and statistic)

    sns.set(style="whitegrid")

    normed_analyses, md = normalize_analyses(analyses)

    df = collate_analysis_results(normed_analyses, field=field, rank=rank)

    stat = []
    if statistic:
        for analysis in normed_analyses:
            try:
                stat.append(alpha_diversity(analysis, statistic, field=field, rank=rank))
            except:
                raise NotImplementedError('{} statistic are currently supported.'.format(statistic))
    elif tax_id:
        stat = df.loc[:, str(tax_id)].values

    md['_data'] = stat

    # numpy types: object, int64, float64, datetime64/timedelta[ns]
    # columns with mixed types are assigned object
    # columns with numbers and NaNs will default to float64
    if 'date' in metadata.split('_'):
        md.loc[:, metadata] = md.loc[:, metadata].apply(pd.to_datetime, utc=True)
        fig, ax = plt.subplots()
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
        ax.plot_date(md[metadata].values, md['_data'].values)
        fig.autofmt_xdate()
    elif md[metadata].dtype == np.object:
        md.fillna(value='N/A', inplace=True)
        sns.boxplot(x=metadata, y='_data', data=md, palette='pastel')
    else:
        sns.lmplot(x=metadata, y='_data', truncate=True, data=md)

    if statistic:
        plt.gca().set_ylabel(statistic)
    elif tax_id:
        plt.gca().set_ylabel(field)

    _, labels = plt.xticks()
    plt.setp(labels, rotation=90)

    if title:
        plt.suptitle(title)

    plt.show()
