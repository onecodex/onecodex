import pandas as pd

# from onecodex.exceptions import OneCodexException
from onecodex.viz.helpers import normalize_analyses, collate_analysis_results


def plot_metadata(analyses, title=None, metadata='created_at', statistic=None,
                  tax_id=None, field='readcount_w_children', rank='species'):
    # metadata -> string (metadata field in default *or* custom metadata) -> x axis?
    # metadata -> (lambda?) # TODO
    # ONE of:
    # tax_id -> plot abundance on y axis
    # statistic -> plot statistic on y axis # TODO
    # statistic (lambda?) # TODO
    import matplotlib.pyplot as plt
    import seaborn as sns
    assert tax_id or statistic
    assert not tax_id and statistic

    sns.set(style="whitegrid")

    normed_analyses, metadata_objs = normalize_analyses(analyses)
    df = collate_analysis_results(normed_analyses, field=field)
    df = df.loc[:, [i[2] == rank for i in df.columns]]
    if tax_id:
        df = df.loc[:, str(tax_id)]

    md = {}
    for idx, analysis in enumerate(normed_analyses):
        if hasattr(metadata_objs[idx], metadata):
            md[analysis.id] = getattr(metadata_objs[idx], metadata)
        elif metadata in metadata_objs[idx].custom:
            md[analysis.id] = metadata_objs[idx].custom[metadata]

    md = pd.DataFrame(md, index=[0])
    md.index = [metadata]

    df = pd.concat([df[:].T, md]).T
    df.rename(columns={df.columns.values[0]: field}, inplace=True)

    sns.boxplot(x=metadata, y=field, data=df, palette='vlag')
    sns.swarmplot(x=metadata, y=field, data=df, size=2, color=".3", linewidth=0)
    plt.show()
