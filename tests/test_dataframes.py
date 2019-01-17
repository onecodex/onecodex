import pytest; pytest.importorskip('pandas')  # noqa
import pandas as pd

from onecodex.analyses import AnalysisMixin
from onecodex.dataframes import ClassificationsDataFrame, ClassificationsSeries, OneCodexAccessor


def test_pandas_subclass():
    inner_df = pd.DataFrame({'datum1': [7, 4, 21], 'datum2': [8, 16, 24]})

    ocx_data = {
        'ocx_rank': 'jedi',
        'ocx_field': 'bantha poodoo',
        'ocx_taxonomy': inner_df.copy(),
        'ocx_metadata': inner_df.copy()
    }

    df = ClassificationsDataFrame({'1279': [1, 2, 3], '1280': [4, 5, 6]}, **ocx_data)

    # we want to be sure that slices of our df are returned as ResultsSeries
    assert type(df['1279']) is ClassificationsSeries

    # we're mostly interested in whether our metadata is transferred between copies. some operations
    # split the df into series and concat back to df, so by doing all this stuff to it we're actually
    # testing several consecutive copy operations using different parts of the pandas API
    new_df = (((df * 10) / 2.2).round(5)) ** 0.5

    # rank is explicitly /not/ passed on, since we don't know what the user has done to the df and
    # we therefore can't trust the rank to be correct
    assert new_df.ocx_rank is None
    assert new_df.ocx_field == 'bantha poodoo'
    assert (new_df.ocx_taxonomy == inner_df).all().all()
    assert (new_df.ocx_metadata == inner_df).all().all()


def test_pandas_extension(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')
    results = samples.to_df()

    # extension should be in ocx namespace of ClassificationsDataFrame
    assert getattr(results, 'ocx', False)
    assert isinstance(results.ocx, OneCodexAccessor)
    assert type(results.ocx).__base__ == AnalysisMixin

    # changes to contents of results df should affect contents of taxonomy df, by keeping only
    # tax_ids in the results df and their parents
    results = samples.to_df(top_n=2)
    assert sorted(results.ocx.taxonomy.index.tolist(), key=int) == \
        ['1', '2', '191', '815', '816', '976', '1224', '28211', '41295',
         '68336', '131567', '171549', '200643', '204441', '1783270']
