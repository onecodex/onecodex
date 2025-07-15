from datetime import datetime, timezone
from typing import Annotated, Union

from pydantic import BeforeValidator, Field, PlainSerializer


def rfc3339_encoder(dt: datetime) -> str:
    """Encode a datetime as an RFC3339 timestamp in UTC (Zulu time).

    If the datetime is timezone-naive, it is assumed to be in UTC.
    If the datetime is in a different timezone, it is converted to UTC.
    """
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    else:
        dt = dt.astimezone(timezone.utc)
    return dt.isoformat()


def rfc3339_decoder(dt: Union[str, datetime]) -> datetime:
    """Decode an RFC3339 timestamp to a datetime."""
    if isinstance(dt, str):
        return datetime.fromisoformat(dt)
    return dt


RFC3339Datetime = Annotated[
    datetime,
    PlainSerializer(rfc3339_encoder, return_type=str),
    BeforeValidator(rfc3339_decoder),
    Field(json_schema_extra={"format": "date-time"}),
]
