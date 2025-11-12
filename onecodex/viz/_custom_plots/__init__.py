from __future__ import annotations

from typing import TYPE_CHECKING, Any
import json

from onecodex.utils import get_requests_session
from onecodex.viz import configure_onecodex_theme
from onecodex.exceptions import OneCodexException
from .collection import SampleCollection, Samples
from .models import PlotParams
from .enums import SuggestionType

if TYPE_CHECKING:
    from pyodide.ffi import JsProxy


def init():
    configure_onecodex_theme()


def plot(params: JsProxy, csrf_token: str) -> dict:
    import js  # available from pyodide

    params = PlotParams.model_validate(_convert_jsnull_to_none(params.to_py()))

    # Cache the *unfiltered* SampleCollection, keyed by (<tag_or_project_uuid>, <metric>).
    # It looks gross but this is how you cache in Pyodide ¯\_(ツ)_/¯
    cache = globals().get("CUSTOM_PLOTS_CACHE", {})

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

    key = (uuid, params.metric)
    if key in cache:
        collection = cache[key]
    else:
        base_url = js.self.location.origin
        url = f"{base_url}/api/v2/custom-plots/sample-data"

        samples = []
        page = 1
        while True:
            session = get_requests_session(headers={"X-CSRFToken": csrf_token})
            resp = session.get(url, params={"type": type_, "uuid": uuid, "page": page})
            resp.raise_for_status()
            samples.extend(Samples(sample) for sample in resp.json())

            num_samples = len(samples)
            if num_samples >= json.loads(resp.headers.get("X-Pagination", "{}")).get("total", 0):
                break
            page += 1

        collection = SampleCollection(samples, metric=params.metric)
        cache[key] = collection

    # Now that we have a SampleCollection, do the actual plotting
    result = collection.plot(params)

    # Cache the results in case the original metric was "auto" and got resolved to a concrete metric
    cache[(uuid, result.metric)] = collection
    globals()["CUSTOM_PLOTS_CACHE"] = cache

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
