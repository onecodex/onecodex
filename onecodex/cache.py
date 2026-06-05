"""Cross-process disk cache for immutable analysis result fetches.

Cache entries are keyed by a hash of (sdk_version, resource_type, id, job_id, method,
args) so that re-analysis under a new job version invalidates cleanly and an SDK
upgrade can't deserialize a stale shape. Only entries for completed, successful
analyses are written.

The on-disk format is zstd-compressed orjson — both deps are already required by the
SDK, the bytes are inspectable with `zstdcat | jq`, and the format survives Python
version changes (unlike pickle).
"""

from __future__ import annotations

import functools
import hashlib
import logging
import os
import secrets
import tempfile
from pathlib import Path
from typing import Any, Optional

import orjson
import zstandard

from onecodex.version import __version__

log = logging.getLogger(__name__)


def _default_cache_dir() -> Path:
    # Default to a per-user dir under the system tempdir so the cache doesn't grow
    # unboundedly on long-lived machines: tempdirs are reclaimed across reboots
    # (launchd on macOS, systemd-tmpfiles on Linux). Callers that want persistence
    # across reboots should pass an explicit path instead.
    return Path(tempfile.gettempdir()) / f"onecodex-results-cache-{os.getuid()}"


class DiskCache:
    """Cross-process content-addressed cache for parsed analysis results.

    Used by the `disk_cached` decorator to memoize the return value of analysis
    `_results` / `_filtered_results` methods. The cache stores the parsed dict
    returned by the method, not the raw HTTP/S3 payload — so it works uniformly
    for analyses that fetch from a pre-signed S3 URL (Classifications) and those
    that fetch from the API (Panels, Mlsts, Workflows, FunctionalProfiles, ...).

    Entries are written atomically (`<key>.tmp.<pid>.<rand>` → `os.replace`), so
    multiple processes sharing a cache dir are safe without locking. Corrupt
    entries (bad zstd / bad json — e.g. from a partial-write race or a tool
    truncating files) are unlinked on read and treated as a miss.

    On-disk format: `zstd(orjson(value))` at `<path>/<sha256>.json.zst`. Both
    deps are already required by the SDK, the format survives Python upgrades
    (unlike pickle), and entries are inspectable with `zstdcat <file> | jq`.

    Keys are produced by `_key(...)` and include the SDK version, the resource
    type and id, an `updated_at` freshness stamp (so a cleared-and-rerun
    analysis invalidates cleanly), the method name, and any args/kwargs.
    """

    def __init__(self, path: Optional[str] = None):
        self.path = Path(path) if path else _default_cache_dir()
        self.path.mkdir(parents=True, exist_ok=True)

    def _key(
        self,
        resource_type: str,
        resource_id: str,
        updated_at: str,
        method: str,
        args: tuple,
        kwargs: dict,
    ) -> str:
        # updated_at is the freshness signal: an analysis can be cleared and re-run
        # under the same uuid, so the uuid alone is not a stable handle. The SDK
        # version is part of the key so a release that changes return shape can't
        # deserialize stale entries.
        payload = orjson.dumps(
            {
                "v": __version__,
                "t": resource_type,
                "id": resource_id,
                "c": updated_at,
                "m": method,
                "a": args,
                "k": kwargs,
            },
            option=orjson.OPT_SORT_KEYS,
        )
        return hashlib.sha256(payload).hexdigest()

    def _entry_path(self, key: str) -> Path:
        return self.path / f"{key}.json.zst"

    def get(self, key: str) -> Optional[Any]:
        entry = self._entry_path(key)
        try:
            blob = entry.read_bytes()
        except FileNotFoundError:
            return None
        try:
            return orjson.loads(zstandard.decompress(blob))
        except (zstandard.ZstdError, orjson.JSONDecodeError, ValueError) as e:
            log.warning("onecodex disk cache: dropping corrupt entry %s (%s)", entry, e)
            try:
                entry.unlink()
            except FileNotFoundError:
                pass
            return None

    def set(self, key: str, value: Any) -> None:
        entry = self._entry_path(key)
        tmp = entry.with_suffix(f".tmp.{os.getpid()}.{secrets.token_hex(4)}")
        try:
            tmp.write_bytes(zstandard.compress(orjson.dumps(value)))
            os.replace(tmp, entry)
        except OSError as e:
            log.warning("onecodex disk cache: failed to write %s (%s)", entry, e)
            try:
                tmp.unlink()
            except FileNotFoundError:
                pass


def disk_cached(method):
    """Cache the return value of an `_AnalysesBase` results method on disk.

    Skipped (read and write) when the owning `Api` has no cache configured or the
    analysis isn't successfully complete. An existing analysis can be cleared and
    re-run under the same uuid, so the key includes `updated_at` to invalidate
    cleanly on re-run.
    """

    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        cache = getattr(self._api, "_cache", None)
        if cache is None or not self.complete or not self.success:
            return method(self, *args, **kwargs)
        key = cache._key(
            resource_type=type(self).__name__,
            resource_id=self.id,
            updated_at=self.updated_at.isoformat(),
            method=method.__name__,
            args=args,
            kwargs=kwargs,
        )
        hit = cache.get(key)
        if hit is not None:
            return hit
        value = method(self, *args, **kwargs)
        cache.set(key, value)
        return value

    return wrapper
