from collections import Counter
from typing import Any
from datetime import datetime, date
from numbers import Number
from functools import partial

from onecodex.lib.enums import BaseEnum


def metadata_record_to_label(metadata: dict, label_by: list[str]) -> str:
    values = []
    for field in label_by:
        value = get_metadata_field_value(metadata, field)
        values.append(normalize_value(value))
    return " / ".join(values)


# TODO is this needed?
def normalize_value(value: Any) -> str:
    try:
        # Metadata stores dates as ISO 8601 strings
        if isinstance(value, str):
            value = datetime.fromisoformat(value).date()
    except ValueError:
        pass  # Do not care

    if isinstance(value, datetime):
        value = value.date()

    if isinstance(value, date):
        return value.strftime("%Y-%m-%d")  # US formating would not work nice with sorting
    elif isinstance(value, BaseEnum):
        return value.value

    return str(value)


def get_metadata_field_value(metadata: dict, field: str) -> Any:
    import pandas as pd

    value = metadata.get(field)
    if pd.isnull(value):
        value = "N/A"
    return value


def deduplicate_labels(labels_by_metadata_id: dict) -> dict:
    num_unique_labels = Counter(labels_by_metadata_id.values())
    new_labels_idx = Counter()

    unique_labels_by_metadata_id = {}
    for metadata_id, label in labels_by_metadata_id.items():
        if num_unique_labels[label] > 1:
            new_labels_idx[label] += 1
            unique_label = f"{label} ({new_labels_idx[label]})"
        else:
            unique_label = label
        unique_labels_by_metadata_id[metadata_id] = unique_label

    return unique_labels_by_metadata_id


def sort_metadata_records(metadata_records: list[dict], sort_by: str) -> list[dict]:
    get_value = partial(get_metadata_field_value, field=sort_by)

    all_keys = {get_value(x) for x in metadata_records}
    if all(isinstance(x, Number) for x in all_keys):
        # TODO does this work?
        return sorted(metadata_records, key=get_value)
    else:
        return sorted(metadata_records, key=lambda x: normalize_value(get_value(x)).lower())
