import pytest

pytest.importorskip("pandas")  # noqa
import pandas as pd

from onecodex.analyses import AnalysisMixin
from onecodex.dataframes import ClassificationsDataFrame, ClassificationsSeries, OneCodexAccessor


def test_pandas_subclass():
    inner_df = pd.DataFrame({"datum1": [7, 4, 21], "datum2": [8, 16, 24]})

    ocx_data = {
        "ocx_rank": "jedi",
        "ocx_metric": "bantha poodoo",
        "ocx_taxonomy": inner_df.copy(),
        "ocx_metadata": inner_df.copy(),
    }

    df = ClassificationsDataFrame({"1279": [1, 2, 3], "1280": [4, 5, 6]}, **ocx_data)

    # we want to be sure that slices of our df are returned as ResultsSeries
    assert type(df["1279"]) is ClassificationsSeries

    # We're mostly interested in whether our metadata is transferred between copies. Here we're just
    # testing that the metadata carries over between slicing. There's another test below that tests
    # other mathematical operations. There is currently a bug upstream in pandas where subclass
    # metadata doesn't carry through certain operations, but does on others.
    new_df = df[["1279"]]

    # rank is explicitly /not/ passed on, since we don't know what the user has done to the df and
    # we therefore can't trust the rank to be correct
    assert new_df.ocx_rank == "jedi"
    assert new_df.ocx_metric == "bantha poodoo"
    assert (new_df.ocx_taxonomy == inner_df).all().all()
    assert (new_df.ocx_metadata == inner_df).all().all()


def test_pandas_extension(samples):
    samples._collate_results(metric="readcount_w_children")
    results = samples.to_df()

    # extension should be in ocx namespace of ClassificationsDataFrame
    assert getattr(results, "ocx", False)
    assert isinstance(results.ocx, OneCodexAccessor)
    assert type(results.ocx).__base__ == AnalysisMixin

    # changes to contents of results df should affect contents of taxonomy df, by keeping only
    # tax_ids in the results df and their parents
    results = samples.to_df(top_n=2, rank="genus")
    assert sorted(results.ocx.taxonomy.index.tolist(), key=int) == [
        "1",
        "2",
        "191",
        "815",
        "816",
        "976",
        "1224",
        "28211",
        "41295",
        "68336",
        "131567",
        "171549",
        "200643",
        "204441",
        "1783270",
    ]


# Current, there's a bug in pandas where these types of operations don't carryover
# the metadata on the subclass.
# https://github.com/pandas-dev/pandas/issues/34177
@pytest.mark.xfail
def test_pandas_subclass_math():
    inner_df = pd.DataFrame({"datum1": [7, 4, 21], "datum2": [8, 16, 24]})

    ocx_data = {
        "ocx_rank": "jedi",
        "ocx_metric": "bantha poodoo",
        "ocx_taxonomy": inner_df.copy(),
        "ocx_metadata": inner_df.copy(),
    }

    df = ClassificationsDataFrame({"1279": [1, 2, 3], "1280": [4, 5, 6]}, **ocx_data)

    # we're mostly interested in whether our metadata is transferred between copies. some operations
    # split the df into series and concat back to df, so by doing all this stuff to it we're actually
    # testing several consecutive copy operations using different parts of the pandas API
    new_df = (((df * 10) / 2.2).round(5)) ** 0.5

    # rank is explicitly /not/ passed on, since we don't know what the user has done to the df and
    # we therefore can't trust the rank to be correct
    assert new_df.ocx_rank == "jedi"
    assert new_df.ocx_metric == "bantha poodoo"
    assert (new_df.ocx_taxonomy == inner_df).all().all()
    assert (new_df.ocx_metadata == inner_df).all().all()
