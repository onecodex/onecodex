import mock
import pytest

from onecodex.exceptions import ValidationError
from onecodex.metadata_upload import (
    validate_appendables,
    validate_tags,
    validate_metadata,
    validate_metadata_against_schema,
    validate_number,
    validate_enum,
    validate_boolean,
    validate_datetime,
    is_blacklisted,
    coerce_custom_value,
)


def schema_rules(ocx, value):
    if value:
        return ocx.Metadata._resource._schema["properties"][value]
    else:
        return ocx.Metadata._resource._schema["properties"]


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
    too_many_characters = ocx.Tags._resource._schema["properties"]["name"]["maxLength"] + 1
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
    schema_props = schema_rules(ocx, None)
    with mock.patch("onecodex.metadata_upload.validate_enum") as mock_enum:
        validate_metadata_against_schema(schema_props, "platform", "Illumina")
        assert mock_enum.call_count == 1

    with mock.patch("onecodex.metadata_upload.validate_number") as mock_number:
        validate_metadata_against_schema(schema_props, "location_lat", "90")
        assert mock_number.call_count == 1

    with mock.patch("onecodex.metadata_upload.validate_boolean") as mock_bool:
        validate_metadata_against_schema(schema_props, "starred", "true")
        assert mock_bool.call_count == 1

    with mock.patch("onecodex.metadata_upload.validate_datetime") as mock_datetime:
        validate_metadata_against_schema(schema_props, "date_collected", "2018,05,08,12,12")
        assert mock_datetime.call_count == 1
    string_response = validate_metadata_against_schema(schema_props, "name", "foo")
    assert string_response == "foo"


def test_validate_number_valid(ocx):
    schema_props = schema_rules(ocx, "location_lat")
    valid_number = validate_number("50", schema_props)
    assert valid_number == 50


def test_validate_number_invalid_large(ocx):
    schema_props = schema_rules(ocx, "location_lat")
    with pytest.raises(ValidationError) as exception_info:
        validate_number("200", schema_props)
    assert "200 must be smaller than the maximum value: 90.0" in str(exception_info)


def test_validate_number_invalid_small(ocx):
    schema_props = schema_rules(ocx, "location_lat")
    with pytest.raises(ValidationError) as exception_info:
        validate_number("-200", schema_props)
    assert "-200 must be larger than the minimum value: -90.0" in str(exception_info)


def test_validate_enum_valid(ocx):
    schema_props = schema_rules(ocx, "platform")
    validate_enum("Illumina HiSeq", schema_props)


def test_validate_enum_invalid(ocx):
    schema_props = schema_rules(ocx, "platform")
    with pytest.raises(ValidationError) as exception_info:
        validate_enum("Foo", schema_props)
    assert "Foo is not a valid value for this key." in str(exception_info)


def test_validate_boolean_truthy():
    validate_boolean("true")
    validate_boolean("TRUE")
    validate_boolean("T")
    validate_boolean("1")
    validate_boolean("Y")
    validate_boolean("YES")


def test_validate_boolean_falsy():
    validate_boolean("false")
    validate_boolean("FALSE")
    validate_boolean("F")
    validate_boolean("0")
    validate_boolean("N")
    validate_boolean("No")


def test_validate_boolean_invalid():
    with pytest.raises(ValidationError) as exception_info:
        validate_boolean("FOO")
    assert 'FOO must be either "true" or "false"' in str(exception_info)


def test_validate_datetime_valid():
    validate_datetime("2018-05-15T16:21:36+00:00")
    validate_datetime("2018-05-15T16:21:36")


def test_validate_datetime_invalid():
    with pytest.raises(ValidationError) as exception_info:
        validate_datetime("2018, 05, 15, 16, 21, 36")
    assert (
        '"2018, 05, 15, 16, 21, 36" must be formatted in iso8601 compliant date format. Example: "2018-05-15T16:21:36+00:00"'
        in str(exception_info)
    )


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


def test_pandas_and_numpy_type_coercions(ocx, upload_mocks):
    # See https://github.com/onecodex/onecodex/issues/232 and
    # https://github.com/pandas-dev/pandas/issues/25969
    # We can remove this test once the Pandas bug above is fixed
    # upstream and we require >= that version of pandas.
    pytest.importorskip("pandas")
    import numpy as np
    import pandas as pd

    series = pd.Series({"a": np.int64(64), "b": 10, "c": "ABC"})
    metadata = {
        "platform": "Illumina NovaSeq 6000",
        "date_collected": "2019-04-14T00:51:54.832048+00:00",
        "external_sample_id": "my-lims-ID-or-similar",
        "custom": series.to_dict(),
    }
    # This works
    ocx.Samples._resource.init_multipart_upload(
        filename="SRR2352185.fastq.gz", size=181687821, metadata=metadata
    )

    # Now try to serialize something really not supported, so we can get the Exception
    with pytest.raises(TypeError):
        metadata["custom"]["bad_field"] = ocx.Samples  # not JSON serializable
        ocx.Samples._resource.init_multipart_upload(
            filename="SRR2352185.fastq.gz", size=181687821, metadata=metadata
        )
