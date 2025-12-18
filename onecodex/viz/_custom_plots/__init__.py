from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Tuple
from urllib.parse import urlencode
import orjson

from onecodex.viz import configure_onecodex_theme
from onecodex.exceptions import OneCodexException
from .collection import SampleCollection, Samples
from .enums import PlotType, SamplesFilter, SuggestionType
from .models import PlotParams
from .utils import AsyncRateLimiter

if TYPE_CHECKING:
    from pyodide.ffi import JsProxy
    from pyodide.http import FetchResponse

CUSTOM_PLOTS_CACHE = {}
RESULTS_BATCH_SIZE = 15


def init():
    configure_onecodex_theme()


async def plot(
    params: JsProxy,
    csrf_token: str,
    progress_callback: Callable[[str, float], None] = lambda msg, pct: None,
) -> dict:
    params = PlotParams.model_validate(_convert_jsnull_to_none(params.to_py()))

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

    filter_ = (
        SamplesFilter.WithFunctionalResults
        if params.plot_type == PlotType.Functional
        else SamplesFilter.WithClassifications
    )

    # Cache the SampleCollection, keyed by (<tag_or_project_uuid>, <filter>, <metric>)
    key = (uuid, filter_, params.metric)
    if key in CUSTOM_PLOTS_CACHE:
        collection = CUSTOM_PLOTS_CACHE[key]
    else:
        samples = await _fetch_samples(
            type_=type_,
            uuid=uuid,
            filter_=filter_,
            csrf_token=csrf_token,
            progress_callback=progress_callback,
        )
        collection = SampleCollection(samples, metric=params.metric)
        CUSTOM_PLOTS_CACHE[key] = collection

    # Now that we have a SampleCollection, do the actual plotting
    result = collection.plot(params)

    # Cache the results in case the original metric was "auto" and got resolved to a concrete metric
    CUSTOM_PLOTS_CACHE[(uuid, filter_, result.params.metric)] = collection

    return result.to_dict()


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


async def _fetch_ndjson(url: str, headers: dict | None = None) -> list[dict[str, Any]]:
    from pyodide.http import HttpStatusError

    resp = await _fetch_with_retries(
        url=url, headers=headers, retries=5, status_forcelist=(202, 429, 502, 503)
    )

    if resp.status == 202:
        # Means data is not ready
        raise HttpStatusError(status=202, status_text="Analysis results not ready", url=url)
    resp.raise_for_status()

    text = (
        await resp.string()
    )  # NDJSON (already decompressed if gzip since browsers do that automatically)
    return [orjson.loads(line) for line in text.splitlines() if line.strip()]


async def _fetch_batch_sample_results(base_url: str, data: list[dict], headers: dict | None = None):
    """
    Download a batch of analysis results.

    Downloads either classification or functional results for a batch of samples depending on
    which field is populated.
    """
    classification_uuids = []
    functional_uuids = []
    # Mapping analysis uuids to position in data list
    uuid_pos = {}

    for idx, sample_data in enumerate(data):
        if sample_data.get("primary_classification"):
            uuid = sample_data["primary_classification"].get("uuid")
            classification_uuids.append(uuid)
            uuid_pos[uuid] = idx
        elif sample_data.get("functional_profile"):
            uuid = sample_data["functional_profile"].get("uuid")
            functional_uuids.append(uuid)
            uuid_pos[uuid] = idx

    if classification_uuids:
        from .utils import format_classification_results

        url = f"{base_url}/api/v2/custom-plots/classification-results?uuids={','.join(classification_uuids)}"
        # Results are ordered according to order in the uuid list
        results = await _fetch_ndjson(url, headers=headers)
        for uuid, result in zip(classification_uuids, results):
            data[uuid_pos[uuid]]["primary_classification"]["api_results"] = (
                format_classification_results(result)
            )

    if functional_uuids:
        url = (
            f"{base_url}/api/v2/custom-plots/functional-results?uuids={','.join(functional_uuids)}"
        )
        results = await _fetch_ndjson(url, headers=headers)
        for uuid, result in zip(functional_uuids, results):
            data[uuid_pos[uuid]]["functional_profile"]["results"] = result


async def _fetch_results_for_samples(base_url: str, data: list[dict], headers: dict | None = None):
    # 5 batch actions per second
    limiter = AsyncRateLimiter(max_actions=5, period=1.0)

    # Rate limited batch streaming sample results
    for i in range(0, len(data), RESULTS_BATCH_SIZE):
        batch = data[i : i + RESULTS_BATCH_SIZE]
        await limiter.acquire()
        await _fetch_batch_sample_results(base_url, batch, headers)


async def _fetch_samples(
    *,
    type_: SuggestionType,
    uuid: str,
    filter_: SamplesFilter,
    csrf_token: str,
    progress_callback: Callable[[str, float], None] = lambda msg, pct: None,
) -> list[Samples]:
    import js  # available from pyodide

    base_url = js.self.location.origin
    url = f"{base_url}/api/v2/custom-plots/sample-data"
    headers = {
        "X-CSRFToken": csrf_token,
        "Accept": "application/json",
    }

    samples = []
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
        sample_data = await resp.json()

        # Fetch results data
        await _fetch_results_for_samples(base_url=base_url, data=sample_data, headers=headers)

        # Convert to ocx sample collection
        samples.extend(Samples(sample) for sample in sample_data)

        pagination = orjson.loads(resp.headers.get("x-pagination", "{}"))
        total = int(pagination.get("total", 0))
        next_page = int(pagination.get("next_page", 0))
        progress_callback("Loading samples", len(samples) / (total or 1))

    return samples


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
        resp = await asyncio.wait_for(pyfetch(url, method=method, headers=headers), timeout=timeout)
        if attempt < retries and resp.status in status_forcelist:
            delay = backoff_factor * (2**attempt)  # exponential backoff
            await asyncio.sleep(delay)
        else:
            return resp


__all__ = ["init", "plot"]
