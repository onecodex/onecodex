"""Tests for the opt-in disk cache (onecodex.cache, Api(cache_results=...))."""

from unittest import mock

import pytest

from onecodex import Api
from onecodex.cache import DiskCache


@pytest.fixture
def ocx_cached(tmp_path):
    return Api(
        api_key="1eab4217d30d42849dbde0cd1bb94e39",
        base_url="http://localhost:3000",
        cache_results=str(tmp_path / "cache"),
    )


def test_cache_results_disabled_by_default(api_data):
    ocx = Api(api_key="x", base_url="http://localhost:3000")
    assert ocx._cache is None


def test_cache_results_path(tmp_path):
    ocx = Api(api_key="x", base_url="http://localhost:3000", cache_results=str(tmp_path / "c"))
    assert ocx._cache is not None
    assert ocx._cache.path == tmp_path / "c"
    assert ocx._cache.path.is_dir()


def test_cache_results_true_uses_tempdir(tmp_path, monkeypatch):
    import tempfile

    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    ocx = Api(api_key="x", base_url="http://localhost:3000", cache_results=True)
    # Per-user namespace inside the tempdir so multi-user systems don't collide.
    assert ocx._cache.path.parent == tmp_path
    assert ocx._cache.path.name.startswith("onecodex-results-cache-")
    assert ocx._cache.path.is_dir()


def test_env_var_enables_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("ONECODEX_DISK_CACHE", str(tmp_path / "envcache"))
    ocx = Api(api_key="x", base_url="http://localhost:3000")
    assert ocx._cache.path == tmp_path / "envcache"


def test_classification_results_cached(ocx_cached, api_data):
    """Second call should hit disk cache and skip the HTTP fetch."""
    cls = ocx_cached.Classifications.get("45a573fb7833449a")
    first = cls.results()

    # New model instance (no lru_cache hit) - second call must come from disk.
    cls2 = ocx_cached.Classifications.get("45a573fb7833449a")
    with mock.patch.object(
        cls2._client, "get", side_effect=AssertionError("should not hit network")
    ):
        # get() above already fetched the model; only the /results endpoint matters here.
        # The s3 session is patched separately by mock_requests; results_uri is None in
        # fixtures, so /results would normally be called. Asserting no .get call works.
        second = cls2.results()
    assert second == first


def test_no_cache_when_incomplete(ocx_cached, api_data):
    cls = ocx_cached.Classifications.get("45a573fb7833449a")
    # Force complete=False; cache must not write or read.
    object.__setattr__(cls, "complete", False)
    cls.results()
    # Cache dir should be empty (no entries written).
    assert not any(ocx_cached._cache.path.iterdir())


def test_updated_at_invalidates_key(tmp_path):
    cache = DiskCache(path=str(tmp_path))
    k1 = cache._key("Classifications", "abc", "2025-01-01T00:00:00", "_results", (), {})
    k2 = cache._key("Classifications", "abc", "2025-02-01T00:00:00", "_results", (), {})
    assert k1 != k2


def test_corrupt_entry_treated_as_miss(tmp_path):
    cache = DiskCache(path=str(tmp_path))
    key = cache._key("X", "id", "t", "m", (), {})
    cache._entry_path(key).write_bytes(b"not valid zstd")
    assert cache.get(key) is None
    # Corrupt file is unlinked.
    assert not cache._entry_path(key).exists()


def test_set_get_roundtrip(tmp_path):
    cache = DiskCache(path=str(tmp_path))
    key = cache._key("X", "id", "t", "m", (), {})
    value = {"table": [{"a": 1}], "n_reads": 42}
    cache.set(key, value)
    assert cache.get(key) == value


def test_args_kwargs_in_key(tmp_path):
    cache = DiskCache(path=str(tmp_path))
    k1 = cache._key("F", "id", "t", "_filtered", (), {"annotation": "go"})
    k2 = cache._key("F", "id", "t", "_filtered", (), {"annotation": "ko"})
    assert k1 != k2
