"""Guard against drift between each model's ``.where()`` signature and the
actual filterable fields on the model.

Each ``.where()`` declares its filterable fields as keyword-only parameters
defaulting to the ``UNSET`` sentinel. A server-side schema change that ports
to the client's pydantic model could silently leave the ``where()`` signature
stale. This test fails loudly if that happens.

A field is considered filterable iff its annotation (after stripping
``Optional`` / ``Annotated``) is one of:
    - a scalar (str, int, float, bool, datetime, Enum)
    - a reference to another ``OneCodexBase`` subclass
    - a list of refs

Dict-valued fields, nested-model fields, and lists of plain strings are
excluded.
"""

from __future__ import annotations

import inspect
import types
import typing
from datetime import datetime
from enum import Enum

import pytest
from pydantic import BaseModel

from onecodex.models import (
    Analyses,
    Assets,
    Classifications,
    Documents,
    Jobs,
    Metadata,
    Projects,
    Samples,
    Tags,
    Users,
)
from onecodex.models.base import OneCodexBase, UNSET


# Control-flag parameter names that appear on ``.where()`` signatures but
# aren't field filters. Filtered out before comparing against ``model_fields``.
CONTROL_KWARGS = {
    "sort",
    "limit",
    "public",
    "organization",
    "filter",
    "tags",
}


REGISTRY = [
    Samples,
    Metadata,
    Tags,
    Users,
    Projects,
    Jobs,
    Documents,
    Assets,
    Analyses,
    Classifications,
]


# Fields present on the model that are intentionally excluded from
# ``where()``'s signature. Either the server doesn't support filtering on
# them, or a control flag already handles them.
OMITTED: dict[type, set[str]] = {
    Samples: {
        "primary_classification",  # server 500s on filter
        "tags",  # resolved via the ``tags=`` control kwarg
    },
    Projects: {
        "permissions",  # list of plain strings; not filterable
        "public",  # handled by ``public=`` control flag
    },
    Jobs: {
        "job_args_schema",  # dict; not filterable
        "public",  # handled by ``public=`` control flag
        # Job-creation-only fields not exposed in the server's response schema:
        "cpu",
        "ram_gb",
        "storage_gb",
        "script",
        "image_uri",
        "job_type",
        "description",
        "inject_bearer_token",
    },
    Documents: {
        "downloaders",  # server 500s on filter
    },
    Metadata: {
        "custom",  # arbitrary dict; filter client-side
    },
    Analyses: {
        "dependencies",  # server 500s on filter
        "job_args",  # dict comparator not supported
        "cost",  # nested-model comparator not supported
    },
    Classifications: {
        "dependencies",
        "job_args",
        "cost",
        # Per-request signed S3 URL — filtering on it is meaningless.
        "results_uri",
    },
}


def _unwrap_optional(tp):
    if typing.get_origin(tp) is typing.Annotated or hasattr(tp, "__metadata__"):
        return _unwrap_optional(typing.get_args(tp)[0])
    origin = typing.get_origin(tp)
    if origin is typing.Union or origin is types.UnionType:
        args = [a for a in typing.get_args(tp) if a is not type(None)]
        if len(args) == 1:
            return _unwrap_optional(args[0])
        return tp
    return tp


def _is_filterable(annotation) -> bool:
    inner = _unwrap_optional(annotation)
    origin = typing.get_origin(inner)

    if inner in (str, int, float, bool, datetime):
        return True
    if isinstance(inner, type) and issubclass(inner, Enum):
        return True
    if isinstance(inner, type) and issubclass(inner, OneCodexBase):
        return True
    if origin is typing.Union or origin is types.UnionType:
        return any(_is_filterable(a) for a in typing.get_args(inner))
    if origin in (list, typing.List):
        (elem,) = typing.get_args(inner) or (object,)
        elem = _unwrap_optional(elem)
        elem_origin = typing.get_origin(elem)
        if elem_origin is typing.Union or elem_origin is types.UnionType:
            return any(
                isinstance(a, type) and issubclass(a, OneCodexBase) for a in typing.get_args(elem)
            )
        return isinstance(elem, type) and issubclass(elem, OneCodexBase)

    return False


def _filterable_fields(model: type[BaseModel]) -> set[str]:
    """Walk ``model_fields`` and return the names of fields that match the
    server's filterable-type criteria. ``field_uri`` is excluded: users
    address it as ``id`` and the client transparently rewrites the kwarg."""
    return {
        name
        for name, info in model.model_fields.items()
        if name != "field_uri" and _is_filterable(info.annotation)
    }


def _declared_filter_kwargs(model: type) -> set[str]:
    """Return the names of keyword-only params on ``model.where()`` that
    default to ``UNSET`` — i.e., the inline filter-field params."""
    sig = inspect.signature(model.where)
    return {
        name
        for name, p in sig.parameters.items()
        if p.kind is inspect.Parameter.KEYWORD_ONLY
        and p.default is UNSET
        and name not in CONTROL_KWARGS
    }


@pytest.mark.parametrize("model", REGISTRY, ids=lambda m: m.__name__)
def test_where_signature_matches_schema(model):
    """Every filterable field on the model must either appear as a keyword-only
    param on ``where()`` OR be explicitly listed in ``OMITTED``."""
    schema_fields = _filterable_fields(model)
    declared = _declared_filter_kwargs(model)
    omitted = OMITTED.get(model, set())

    missing = schema_fields - declared - omitted
    extra = declared - schema_fields

    assert not missing, (
        f"{model.__name__}: filterable fields missing from .where() signature: "
        f"{sorted(missing)}. Either add as a keyword-only param or to OMITTED "
        f"with a reason."
    )
    assert not extra, (
        f"{model.__name__}: .where() declares unknown field(s): {sorted(extra)}. "
        f"Remove from the signature or check the model schema."
    )
