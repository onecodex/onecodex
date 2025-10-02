import pytest

from onecodex.exceptions import ValidationError
from onecodex.metadata_upload import (
    validate_appendables,
    validate_tags,
    validate_metadata,
    validate_metadata_against_schema,
    is_blacklisted,
    coerce_custom_value,
)


def test_validate_appendables(ocx):
    initial_appendables = {"metadata": {"foo": "bar", "starred": "true"}, "tags": ["baz"]}
    final_appendables = dict(initial_appendables)
    final_appendables["valid_metadata"] = {"starred": True}
    final_appendables["valid_metadata"]["custom"] = {"foo": "bar"}
    final_appendables["valid_tags"] = [{"name": "baz"}]
    appendables = validate_appendables(initial_appendables, ocx)
    assert appendables == final_appendables


def test_validate_tags_valid(ocx):
    initial_appendables = {
        "metadata": {},
        "tags": ["baz"],
        "valid_tags": [],
        "valid_metadata": {"custom": {}},
    }
    validate_tags(initial_appendables, ocx)
    assert initial_appendables["valid_tags"] == [{"name": "baz"}]


def test_validate_tags_invalid(ocx):
    too_many_characters = ocx.Tags.model_fields["name"].metadata[0].max_length + 1
    invalid_tag = "a" * too_many_characters
    initial_appendables = {
        "metadata": {},
        "tags": [invalid_tag],
        "valid_tags": [],
        "valid_metadata": {"custom": {}},
    }
    with pytest.raises(ValidationError):
        validate_tags(initial_appendables, ocx)


def test_validate_metadata_valid(ocx):
    initial_appendables = {
        "metadata": {"starred": "true"},
        "tags": [],
        "valid_tags": [],
        "valid_metadata": {"custom": {}},
    }
    final_appendables = dict(initial_appendables)
    final_appendables["valid_metadata"]["starred"] = True
    validate_metadata(initial_appendables, ocx)
    assert initial_appendables == final_appendables


def test_validate_metadata_not_present(ocx):
    initial_appendables = {"tags": [], "valid_tags": [], "valid_metadata": {"custom": {}}}
    final_appendables = dict(initial_appendables)
    validate_metadata(initial_appendables, ocx)
    assert initial_appendables == final_appendables


def test_validate_metadata_blacklisted(ocx):
    initial_appendables = {
        "metadata": {"$uri": "invalid_entry"},
        "tags": [],
        "valid_tags": [],
        "valid_metadata": {"custom": {}},
    }
    with pytest.raises(ValidationError):
        validate_metadata(initial_appendables, ocx)


def test_validate_metadata_custom_valid(ocx):
    initial_appendables = {
        "metadata": {"foo": "bar"},
        "tags": [],
        "valid_tags": [],
        "valid_metadata": {"custom": {}},
    }
    final_appendables = dict(initial_appendables)
    final_appendables["valid_metadata"]["custom"]["foo"] = "bar"
    validate_metadata(initial_appendables, ocx)
    assert initial_appendables == final_appendables


def test_validate_metadata_against_schema(ocx):
    validate_metadata_against_schema("platform", "Illumina")
    validate_metadata_against_schema("location_lat", "90")
    validate_metadata_against_schema("starred", "true")
    with pytest.raises(ValidationError):
        validate_metadata_against_schema("date_collected", "2018,05,08,12,12")

    string_response = validate_metadata_against_schema("name", "foo")
    assert string_response == "foo"


def test_validate_number_valid(ocx):
    valid_number = validate_metadata_against_schema("location_lat", "50")
    assert valid_number == 50


def test_validate_number_invalid_large(ocx):
    with pytest.raises(ValidationError) as exception_info:
        validate_metadata_against_schema("location_lat", "200")
    assert "Input should be less" in str(exception_info)


def test_validate_number_invalid_small(ocx):
    with pytest.raises(ValidationError) as exception_info:
        validate_metadata_against_schema("location_lat", "-200")
    assert "Input should be greater" in str(exception_info)


def test_validate_enum_valid(ocx):
    validate_metadata_against_schema("platform", "Illumina HiSeq")


@pytest.mark.xfail(reason="Enum validation is not yet implemented")
def test_validate_enum_invalid(ocx):
    with pytest.raises(ValidationError) as exception_info:
        validate_metadata_against_schema("platform", "Foo")
    assert "Input should be '454 sequencing'" in str(exception_info)


def test_validate_boolean_truthy():
    assert validate_metadata_against_schema("starred", "true") is True
    assert validate_metadata_against_schema("starred", "TRUE") is True
    assert validate_metadata_against_schema("starred", "T") is True
    assert validate_metadata_against_schema("starred", "1") is True
    assert validate_metadata_against_schema("starred", "Y") is True
    assert validate_metadata_against_schema("starred", "YES") is True


def test_validate_boolean_falsy():
    validate_metadata_against_schema("starred", "false") is False
    validate_metadata_against_schema("starred", "FALSE") is False
    validate_metadata_against_schema("starred", "F") is False
    validate_metadata_against_schema("starred", "0") is False
    validate_metadata_against_schema("starred", "N") is False
    validate_metadata_against_schema("starred", "No") is False


def test_validate_boolean_invalid():
    with pytest.raises(ValidationError) as exception_info:
        validate_metadata_against_schema("starred", "FOO")
    assert "Input should be a valid boolean, unable to interpret input" in str(exception_info)


def test_validate_datetime_valid():
    validate_metadata_against_schema("date_collected", "2018-05-15T16:21:36+00:00")
    validate_metadata_against_schema("date_sequenced", "2018-05-15T16:21:36")


def test_validate_datetime_invalid():
    with pytest.raises(ValidationError) as exception_info:
        validate_metadata_against_schema("date_collected", "2018, 05, 15, 16, 21, 36")
    assert "Invalid isoformat string:" in str(exception_info)


def test_is_blacklisted_valid():
    blacklisted_resp = is_blacklisted("platform")
    assert blacklisted_resp is False


def test_is_blacklisted_invalid():
    blacklisted_resp = is_blacklisted("$uri")
    assert blacklisted_resp is True


def test_coerce_custom_value_float():
    custom_value = coerce_custom_value("42")
    assert custom_value == 42.0
    assert type(custom_value) is float


def test_coerce_custom_value_truthy():
    custom_value = coerce_custom_value("True")
    assert custom_value is True


def test_coerce_custom_value_falsy():
    custom_value = coerce_custom_value("False")
    assert custom_value is False


def test_coerce_custom_value_string():
    custom_value = coerce_custom_value("Foo")
    assert custom_value == "Foo"
