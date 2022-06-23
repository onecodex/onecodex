# import json

import pandas as pd

# import pytest

# pytest.importorskip("pandas")  # noqa ||| safay: not sure why this is here
# import warnings

from onecodex.models import SampleCollection
from onecodex.models.experimental import FunctionalProfiles


def test_query_for_functional_analysis(ocx_experimental, api_data):
    sample_id = "73b8349a30b04957"
    profile = ocx_experimental.FunctionalProfiles.where(sample=sample_id)
    assert isinstance(profile[0], FunctionalProfiles)


def test_functional_profiles_table(ocx_experimental, api_data):
    func_profile = ocx_experimental.FunctionalProfiles.get("31ddae978aff475f")
    df = func_profile.table()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 992
    assert set(df.columns) == {'group_name', 'id', 'metric', 'name', 'taxa_stratified', 'taxon_id', 'taxon_name', 'value'}
    assert len(df['taxon_name'].unique()) == 47

    eggnog_df = func_profile.table(annotation="eggnog", taxa_stratified=False)
    assert set(eggnog_df['group_name']) == {'eggnog'}
    assert list(eggnog_df['taxon_name'].unique()) == [None]


def test_functional_profiles_fetch(ocx_experimental, api_data):
    sample_ids = ["543c9c046e3e4e09", "66c1531cb0b244f6", "37e5151e7bcb4f87"]
    samples = [ocx_experimental.Samples.get(sample_id) for sample_id in sample_ids]
    sc = SampleCollection(samples)
    # SampleCollection._functional_profiles_fetch() populates the .functional_profiles attribute cache
    functional_profiles = sc._functional_profiles
    for profile in functional_profiles:
        assert isinstance(profile, FunctionalProfiles)


def test_collate_functional_results(ocx_experimental, api_data):
    sample_ids = ["543c9c046e3e4e09", "66c1531cb0b244f6", "37e5151e7bcb4f87"]
    samples = [ocx_experimental.Samples.get(sample_id) for sample_id in sample_ids]
    sc = SampleCollection(samples)
    df = sc._collate_functional_results(annotation="go", metric="rpk")
    assert isinstance(df, pd.DataFrame)
    print(df.head())
    # # TODO: assert dataframe has expected data
    # TODO: assert some failure paths
