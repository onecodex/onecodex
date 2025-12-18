import pytest

pytest.importorskip("pandas")  # noqa
import pandas as pd

from onecodex.dataframes import ClassificationsDataFrame, ClassificationsSeries


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
