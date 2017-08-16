import pandas as pd

from onecodex.exceptions import OneCodexException
from onecodex.viz.helpers import normalize_analyses, collate_analysis_results
from onecodex.lib.distance_metrics import simpson, chao1


def plot_metadata(analyses, title=None, metadata='created_at', statistic=None,
                  tax_id=None, field='readcount_w_children', rank='species'):
    # metadata -> string (metadata field in default *or* custom metadata) -> x axis?
    # metadata -> (lambda?) # TODO
    # ONE of:
    # tax_id -> plot abundance on y axis
    # statistic -> plot statistic on y axis
    # statistic (lambda?) # TODO
    import matplotlib.pyplot as plt
    import seaborn as sns

    assert tax_id or statistic
    assert not (tax_id and statistic)

    if statistic:
        if statistic == 'chao1':
            calculate = chao1
        elif statistic == 'simpson':
            calculate = simpson

    sns.set(style="whitegrid")
    f, ax = plt.subplots()

    normed_analyses, metadata_objs = normalize_analyses(analyses)

    # FIXME: delete after testing
    for idx, m in enumerate(metadata_objs):
        m.custom['height'] = idx * 15

    df = collate_analysis_results(normed_analyses, field=field)
    df = df.loc[:, [i[2] == rank for i in df.columns]]

    md = {}
    stat = {}
    for idx, analysis in enumerate(normed_analyses):
        if hasattr(metadata_objs[idx], metadata):
            md[analysis.id] = getattr(metadata_objs[idx], metadata)
        elif metadata in metadata_objs[idx].custom:
            md[analysis.id] = metadata_objs[idx].custom[metadata]
        else:
            raise OneCodexException('"{}" is not set for all samples'.format(metadata))

        try:
            md[analysis.id] = float(md[analysis.id])
        except ValueError:
            raise OneCodexException('"{}" cannot be parsed as a number'.format(metadata))

        if statistic:
            stat[analysis.id] = calculate(analysis, field=field, rank=rank)

    if tax_id:
        md = pd.DataFrame(md, index=[0])
        md.index = [metadata]
        df = df.loc[:, str(tax_id)]
        df = pd.concat([df[:].T, md]).T
        df.rename(columns={df.columns.values[0]: field}, inplace=True)
        sns.boxplot(x=metadata, y=field, data=df, palette='vlag', ax=ax)
        sns.swarmplot(x=metadata, y=field, data=df, size=2, color=".3", linewidth=0, ax=ax)
    elif statistic:
        df = pd.DataFrame([md, stat], index=[metadata, statistic]).T
        df.T.index = [metadata, statistic]
        sns.boxplot(x=metadata, y=statistic, data=df.reset_index(), palette='vlag', ax=ax)
        sns.swarmplot(x=metadata, y=statistic, data=df.reset_index(),
                      size=2, color=".3", linewidth=0, ax=ax)

    if title:
        f.suptitle(title)
    plt.show()
