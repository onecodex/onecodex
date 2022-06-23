# import json

import pandas as pd

# import pytest

# pytest.importorskip("pandas")  # noqa ||| safay: not sure why this is here
# import warnings

# from onecodex.models.experimental import FunctionalProfiles


def test_functional_profiles(ocx_experimental, api_data):
    # TODO: add the sample to models
    # sample = ocx_experimental.Samples.get("73b8349a30b04957")
    # func_profile = ocx_experimental.Analysis.where(sample=sample)
    func_profile = ocx_experimental.FunctionalProfiles.get("31ddae978aff475f")
    df = func_profile.table()
    assert isinstance(df, pd.DataFrame)
    print(len(df))
    # TODO: assert the dataframe has expected data...
