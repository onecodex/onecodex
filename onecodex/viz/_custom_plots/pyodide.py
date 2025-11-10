from __future__ import annotations

from typing import TYPE_CHECKING, Any
import json

import requests

from onecodex.viz import configure_onecodex_theme
from .collection import SampleCollection, Samples
from .models import PlotParams
from .enums import SuggestionType

if TYPE_CHECKING:
    from pyodide.ffi import JsProxy


def init():
    configure_onecodex_theme()


def plot(params: JsProxy, csrf_token: str, include_exported_chart_data: bool = False) -> dict:
    import js  # available from pyodide

    params = PlotParams.model_validate(_replace_jsnull(params.to_py()))
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
        raise NotImplementedError

    key = (uuid, params.metric)
    if key in cache:
        collection = cache[key]
    else:
        base_url = js.self.location.origin
        url = f"{base_url}/api/v2/custom-plots/sample-data"

        samples = []
        page = 1
        while True:
            resp = requests.get(
                url,
                params={"type": type_, "uuid": uuid, "page": page},
                headers={"X-CSRFToken": csrf_token},
            )
            # TODO handle 429 retries
            resp.raise_for_status()
            samples.extend(Samples(sample) for sample in resp.json())

            num_samples = len(samples)
            if num_samples >= json.loads(resp.headers.get("X-Pagination", "{}")).get("total", 0):
                break
            page += 1

        collection = SampleCollection(samples, metric=params.metric)
        cache[key] = collection

    result = collection.plot(params)
    cache[(uuid, result.metric)] = collection
    globals()["CUSTOM_PLOTS_CACHE"] = cache

    return result.to_dict(params, include_exported_chart_data=include_exported_chart_data)


def _replace_jsnull(obj: Any) -> Any:
    from pyodide.ffi import jsnull

    if obj is jsnull:
        return None
    elif isinstance(obj, list):
        return [_replace_jsnull(x) for x in obj]
    elif isinstance(obj, dict):
        return {k: _replace_jsnull(v) for k, v in obj.items()}
    return obj
