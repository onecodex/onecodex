from functools import partial

import mock
import pytest

from click import BadParameter

from onecodex.api import Api
from onecodex.utils import (
    snake_case,
    check_for_allowed_file,
    valid_api_key,
    has_missing_values,
    init_sentry,
)


def test_check_allowed_file():
    # bad ones
    with pytest.raises(SystemExit):
        check_for_allowed_file("file.bam")
        check_for_allowed_file("file")

    # good ones
    check_for_allowed_file("file.fastq")
    check_for_allowed_file("file.fastq.gz")


def test_is_valid_api_key():
    empty_key = ""
    short_key = "123"
    long_key = "123abc123abc123abc123abc123abc123abc123abc123abc123abc"
    good_key = "123abc123abc123abc123abc123abc32"

    # its a click callback so it expects some other params
    valid_api_key_partial = partial(valid_api_key, None, None)

    for key in [empty_key, short_key, long_key]:
        with pytest.raises(BadParameter):
            valid_api_key_partial(key)

    assert good_key == valid_api_key_partial(good_key)


@pytest.mark.parametrize(
    "resource,uris",
    [
        ("Samples", []),
        ("Samples", ["761bc54b97f64980"]),
        ("Analyses", []),
        ("Analyses", ["45a573fb7833449a"]),
        ("Markerpanels", []),
    ],
)
def test_fetcher(ocx, api_data, resource, uris):
    if len(uris) == 0:
        pass
    else:
        for uri in uris:
            resource_class = getattr(ocx, resource)
            instance = resource_class.get(uri)
            assert instance is not None


def test_snake_case():
    test_cases = ["SnakeCase", "snakeCase", "SNAKE_CASE"]
    for test_case in test_cases:
        assert snake_case(test_case) == "snake_case"


def test_custom_ca_bundle(runner, api_data):
    """Tests that we're properly merging settings into our prepared requests."""
    with mock.patch("requests.Session.merge_environment_settings") as merge_env:
        ocx = Api(base_url="http://localhost:3000", cache_schema=True)
        classifications = ocx.Classifications.all()
        assert merge_env.call_count >= 1
        assert len(classifications) >= 1


def test_has_missing_values():
    pytest.importorskip("numpy")
    pytest.importorskip("pandas")

    import numpy as np
    import pandas as pd

    assert has_missing_values(pd.Series([1, np.nan, 2]))
    assert has_missing_values(pd.Series([np.nan, np.nan]))
    assert not has_missing_values(pd.Series([1, 2, 3]))

    assert has_missing_values(pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "b", None]}))
    assert not has_missing_values(pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "b", "c"]}))


@pytest.mark.parametrize(
    "ONE_CODEX_NO_TELEMETRY,ONE_CODEX_SENTRY_DSN,call_count,dsn_contains",
    [
        ("1", None, 0, ""),
        (None, None, 1, "sentry.io"),
        (None, "SomeDSN", 1, "SomeDSN"),
    ],
)
def test_init_sentry(
    monkeypatch, ONE_CODEX_NO_TELEMETRY, ONE_CODEX_SENTRY_DSN, call_count, dsn_contains
):
    if ONE_CODEX_NO_TELEMETRY:
        monkeypatch.setenv("ONE_CODEX_NO_TELEMETRY", ONE_CODEX_NO_TELEMETRY)
    if ONE_CODEX_SENTRY_DSN:
        monkeypatch.setenv("ONE_CODEX_SENTRY_DSN", ONE_CODEX_SENTRY_DSN)

    with (
        mock.patch("onecodex.utils._setup_sentry_for_ipython") as _,
        mock.patch("sentry_sdk.init") as mocked_sentry_init,
    ):
        init_sentry()
        assert mocked_sentry_init.call_count == call_count
        if call_count:
            assert dsn_contains in mocked_sentry_init.call_args.kwargs["dsn"]
