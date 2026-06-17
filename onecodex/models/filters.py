"""Operator-shape TypedDicts for ``where()`` filter values.

Each ``.where()`` kwarg accepts either a literal value or a Mongo-style
operator dict. These TypedDicts describe the operator-dict shapes — they're
used as the value type in per-model ``where()`` signatures, and can also be
used by callers who want to type a filter expression they're building::

    from onecodex.models.filters import NumFilter

    big: NumFilter = {"$gte": 1_000_000}
    ocx.Samples.where(size=big)

All TypedDicts use ``total=False`` — every operator key is optional.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TypedDict

from onecodex.models.base import OneCodexBase

# Reference fields accept a model instance or an id string. Single-ref slots
# additionally accept ``None`` (for IS NULL filters).
_Ref = OneCodexBase | str


StrFilter = TypedDict(
    "StrFilter",
    {
        "$eq": str | None,
        "$ne": str | None,
        "$in": Sequence[str],
        "$contains": str,
        "$icontains": str,
        "$startswith": str,
        "$istartswith": str,
        "$endswith": str,
        "$iendswith": str,
    },
    total=False,
)

NumFilter = TypedDict(
    "NumFilter",
    {
        "$eq": int | float | None,
        "$ne": int | float | None,
        "$lt": int | float,
        "$lte": int | float,
        "$gt": int | float,
        "$gte": int | float,
        "$between": Sequence[int | float],
        "$in": Sequence[int | float],
    },
    total=False,
)

# Datetimes go over the wire as RFC3339 strings — ``datetime`` objects aren't
# JSON-serializable by ``json.dumps``. Use ``datetime.isoformat()`` on the
# call site.
DatetimeFilter = TypedDict(
    "DatetimeFilter",
    {
        "$eq": str | None,
        "$ne": str | None,
        "$lt": str,
        "$lte": str,
        "$gt": str,
        "$gte": str,
        "$between": Sequence[str],
        "$in": Sequence[str],
    },
    total=False,
)

BoolFilter = TypedDict(
    "BoolFilter",
    {
        "$eq": bool | None,
        "$ne": bool | None,
    },
    total=False,
)

# For fields where the API only supports equality (no substring, no membership).
# Used on fields like ``Users.email`` where the route's filter validator
# rejects everything but ``$eq`` / ``$ne``.
EqStrFilter = TypedDict(
    "EqStrFilter",
    {
        "$eq": str | None,
        "$ne": str | None,
    },
    total=False,
)

RefFilter = TypedDict(
    "RefFilter",
    {
        "$eq": _Ref | None,
        "$ne": _Ref | None,
        "$in": Sequence[_Ref],
    },
    total=False,
)

ListRefFilter = TypedDict(
    "ListRefFilter",
    {
        "$containsall": Sequence[_Ref],
        "$containsany": Sequence[_Ref],
    },
    total=False,
)


__all__ = [
    "StrFilter",
    "NumFilter",
    "DatetimeFilter",
    "BoolFilter",
    "EqStrFilter",
    "RefFilter",
    "ListRefFilter",
]
