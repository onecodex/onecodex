from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable
import json

from onecodex.utils import get_requests_session
from onecodex.viz import configure_onecodex_theme
from onecodex.exceptions import OneCodexException
from .collection import SampleCollection, Samples
from .enums import SuggestionType
from .models import PlotParams

if TYPE_CHECKING:
    from pyodide.ffi import JsProxy

CUSTOM_PLOTS_CACHE = {}


def init():
    configure_onecodex_theme()


def plot(
    params: JsProxy,
    csrf_token: str,
    progress_callback: Callable[[str, float], None] = lambda msg, pct: None,
) -> dict:
    import js  # available from pyodide

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

    # Cache the *unfiltered* SampleCollection, keyed by (<tag_or_project_uuid>, <metric>)
    key = (uuid, params.metric)
    if key in CUSTOM_PLOTS_CACHE:
        collection = CUSTOM_PLOTS_CACHE[key]
    else:
        base_url = js.self.location.origin
        url = f"{base_url}/api/v2/custom-plots/sample-data"

        samples = []
        next_page = 1
        progress_callback("Loading samples", 0.0)
        while next_page:
            session = get_requests_session(headers={"X-CSRFToken": csrf_token})
            resp = session.get(url, params={"type": type_, "uuid": uuid, "page": next_page})
            resp.raise_for_status()
            samples.extend(Samples(sample) for sample in resp.json())

            pagination = json.loads(resp.headers.get("X-Pagination", "{}"))
            total = int(pagination.get("total", 0))
            next_page = int(pagination.get("next_page", 0))
            progress_callback("Loading samples", len(samples) / (total or 1))

        collection = SampleCollection(samples, metric=params.metric)
        CUSTOM_PLOTS_CACHE[key] = collection

    # Now that we have a SampleCollection, do the actual plotting
    result = collection.plot(params)

    # Cache the results in case the original metric was "auto" and got resolved to a concrete metric
    CUSTOM_PLOTS_CACHE[(uuid, result.params.metric)] = collection

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


__all__ = ["init", "plot"]
