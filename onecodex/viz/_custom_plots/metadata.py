from collections import Counter
from typing import Any


def metadata_record_to_label(metadata: dict, label_by: list[str]) -> str:
    values = []
    for field in label_by:
        values.append(_get_metadata_field_value(metadata, field))
    return " / ".join(values)


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
    return sorted(
        metadata_records, key=lambda record: _get_metadata_field_value(record, sort_by).lower()
    )


def _get_metadata_field_value(metadata: dict, field: str) -> Any:
    import pandas as pd

    value = metadata.get(field)
    if pd.isnull(value):
        value = "N/A"
    return value
