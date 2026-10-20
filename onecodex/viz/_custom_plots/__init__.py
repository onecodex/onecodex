from __future__ import annotations

import asyncio
from typing import TYPE_CHECKING, Any, Callable, Tuple
from urllib.parse import urlencode
import orjson

from onecodex.viz import configure_onecodex_theme
from onecodex.exceptions import ConnectivityError, OneCodexException
from .collection import (
    Classifications,
    FunctionalProfiles,
    SampleCollection,
    samples_from_sample_data,
)
from .enums import PlotType, SamplesFilter, SuggestionType
from .models import BaseParams, PlotParams, StatsParams

if TYPE_CHECKING:
    from pyodide.ffi import JsProxy
    from pyodide.http import FetchResponse

CUSTOM_PLOTS_CACHE = {}

# Keys in the browser's Cache Storage for analysis results downloads
RESULTS_CACHE_PREFIX = "onecodex-custom-plots-results"
RESULTS_CACHE_NAME = f"{RESULTS_CACHE_PREFIX}-v1"
RESULTS_CONCURRENCY = 8


def init():
    configure_onecodex_theme()


async def _get_collection(
    params: BaseParams,
    filter_: SamplesFilter,
    csrf_token: str,
    progress_callback: Callable[[str, float], None],
) -> SampleCollection:
    uuid = None
    type_ = None
    if params.tag:
        uuid = params.tag
        type_ = SuggestionType.Tag
    elif params.project:
        uuid = params.project
        type_ = SuggestionType.Project
    else:
        raise OneCodexException("Neither a tag nor project UUID was provided.")

    # Cache the SampleCollection, keyed by (<tag_or_project_uuid>, <filter>)
    key = (uuid, filter_)
    if key in CUSTOM_PLOTS_CACHE:
        return CUSTOM_PLOTS_CACHE[key]

    collection = await _fetch_collection(
        type_=type_,
        uuid=uuid,
        filter_=filter_,
        csrf_token=csrf_token,
        progress_callback=progress_callback,
    )
    CUSTOM_PLOTS_CACHE[key] = collection
    return collection


async def plot(
    params: JsProxy,
    csrf_token: str,
    progress_callback: Callable[[str, float], None] = lambda msg, pct: None,
) -> dict:
    params = PlotParams.model_validate(_convert_jsnull_to_none(params.to_py()))

    filter_ = (
        SamplesFilter.WithFunctionalResults
        if params.plot_type == PlotType.Functional
        else SamplesFilter.WithClassifications
    )
    collection = await _get_collection(params, filter_, csrf_token, progress_callback)

    return collection.plot(params).to_dict()


async def stats(
    params: JsProxy,
    csrf_token: str,
    progress_callback: Callable[[str, float], None] = lambda msg, pct: None,
) -> dict:
    params = StatsParams.model_validate(_convert_jsnull_to_none(params.to_py()))

    collection = await _get_collection(
        params, SamplesFilter.WithClassifications, csrf_token, progress_callback
    )

    return collection.stats(params).to_dict()


def _convert_jsnull_to_none(obj: Any) -> Any:
    """Convert `jsnull` to `None`.

    Pyodide converts `null` to `jsnull`, and `undefined` to `None`. Convert `jsnull` to `None` too.

    """
    from pyodide.ffi import jsnull

    if obj is jsnull:
        return None
    elif isinstance(obj, list):
        return [_convert_jsnull_to_none(x) for x in obj]
    elif isinstance(obj, dict):
        return {k: _convert_jsnull_to_none(v) for k, v in obj.items()}
    return obj


async def _fetch_collection(
    *,
    type_: SuggestionType,
    uuid: str,
    filter_: SamplesFilter,
    csrf_token: str,
    progress_callback: Callable[[str, float], None] = lambda msg, pct: None,
) -> SampleCollection:
    import js  # available from pyodide

    base_url = js.self.location.origin
    url = f"{base_url}/api/frontend/custom-plots/sample-data"
    headers = {
        "X-CSRFToken": csrf_token,
        "Accept": "application/json",
    }
    results_cache = await _open_results_cache()

    samples = []
    functional_profiles = {}
    next_page = 1
    progress_callback("Loading samples", 0.0)
    while next_page:
        params = urlencode(
            {
                "type": type_,
                "uuid": uuid,
                "filter": filter_,
                "page": next_page,
            }
        )
        full_url = f"{url}?{params}"

        # Fetch sample metadata
        resp = await _fetch_with_retries(url=full_url, headers=headers)
        resp.raise_for_status()
        page_samples, page_functional_profiles = samples_from_sample_data(await resp.json())

        # Fetch results data
        analyses = [s.primary_classification for s in page_samples if s.primary_classification]
        analyses += page_functional_profiles.values()
        await _load_results(analyses, results_cache=results_cache, base_url=base_url)

        samples.extend(page_samples)
        functional_profiles.update(page_functional_profiles)

        pagination = orjson.loads(resp.headers.get("x-pagination", "{}"))
        total = int(pagination.get("total", 0))
        next_page = int(pagination.get("next_page", 0))
        progress_callback("Loading samples", len(samples) / (total or 1))

    return SampleCollection(samples, functional_profiles=functional_profiles)


async def _open_results_cache() -> JsProxy | None:
    try:
        from js import caches

        # Clean up any stale versions
        for name in await caches.keys():
            if name.startswith(RESULTS_CACHE_PREFIX) and name != RESULTS_CACHE_NAME:
                await caches.delete(name)
        return await caches.open(RESULTS_CACHE_NAME)
    except Exception:
        # Cache unavailable, ignoring
        return None


async def _load_results(
    analyses: list[Classifications | FunctionalProfiles],
    *,
    results_cache: JsProxy | None,
    base_url: str,
):
    semaphore = asyncio.Semaphore(RESULTS_CONCURRENCY)

    async def load(analysis: Classifications | FunctionalProfiles):
        async with semaphore:
            body = await _read_results_file(
                analysis, results_cache=results_cache, base_url=base_url
            )
            analysis._loaded_results = orjson.loads(body)

    # Download in parallel with a cap (semaphore)
    await asyncio.gather(*(load(analysis) for analysis in analyses))


async def _read_results_file(
    analysis: Classifications | FunctionalProfiles,
    *,
    results_cache: JsProxy | None,
    base_url: str,
) -> str | bytes:
    """Return the JSON of an analysis' results file, from the results cache if possible."""
    cache_key = f"{base_url}/custom-plots/results/{analysis.id}"

    if results_cache is not None:
        try:
            cached = await results_cache.match(cache_key)
            if cached is not None:
                return await cached.text()
        except Exception:
            # Ignore, just fall back to downloading
            pass

    try:
        resp = await _fetch_with_retries(url=analysis.results_uri)
        resp.raise_for_status()
        # Results should be transparently decompressed
        body = await resp.bytes()
    except AttributeError:
        # Converting pyfetch exceptions to Python exceptions sometimes fails with AttributeError
        # Raising a ConnectivityError that will get handled gracefully
        raise ConnectivityError("Cannot load samples, please try again later")

    if results_cache is not None:
        try:
            from js import Response
            from pyodide.ffi import to_js

            await results_cache.put(cache_key, Response.new(to_js(body)))
        except Exception:
            # Cache isn't critical, ignoring all errors
            pass

    return body


async def _fetch_with_retries(
    *,
    url: str,
    method: str = "GET",
    headers: dict | None = None,
    timeout: int = 30,
    retries: int = 3,
    backoff_factor: float = 4.0,
    status_forcelist: Tuple[int, ...] = (429, 502, 503),
) -> FetchResponse:
    import asyncio
    from pyodide.http import pyfetch

    # We're not using `requests` because it doesn't always work reliably in Pyodide, e.g. with
    # retries or large response payloads. It also generates console errors about setting request
    # headers that are blocked by the browser.
    for attempt in range(retries + 1):
        try:
            resp = await asyncio.wait_for(
                pyfetch(url, method=method, headers=headers), timeout=timeout
            )
            if attempt < retries and resp.status in status_forcelist:
                delay = backoff_factor * (2**attempt)  # exponential backoff
                await asyncio.sleep(delay)
            else:
                return resp
        except TimeoutError:
            if attempt < retries:
                delay = backoff_factor * (2**attempt)  # exponential backoff
                await asyncio.sleep(delay)
            else:
                raise


__all__ = ["init", "plot", "stats"]
