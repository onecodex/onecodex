from onecodex.exceptions import ValidationError, ValidationWarning
from tests.conftest import ocx
import pytest
import mock
from onecodex.metadata_upload import validate_appendables, validate_tags, validate_metadata, validate_metadata_against_schema

def test_validate_appendables():
    initial_appendables = { 'metadata': {'foo': 'bar', 'starred': 'true'}, 'tags': ['baz']}
    final_appendables = dict(initial_appendables)
    final_appendables['custom_metadata'] = {'foo': 'bar'}
    final_appendables['valid_metadata'] = {'starred': True}
    final_appendables['valid_tags'] = ['baz']
    appendables = validate_appendables(initial_appendables, ocx())
    assert appendables == final_appendables


def test_validate_tags_valid():
    initial_appendables = { 'metadata': {}, 'tags': ['baz'], 'valid_tags': [], 'custom_metadata': {}, 'valid_metadata': {}}
    intial_valid_tag_length = len(initial_appendables['valid_tags'])
    validate_tags(initial_appendables, ocx())
    assert len(initial_appendables['valid_tags']) != intial_valid_tag_length


def test_validate_tags_invalid():
    too_many_characters = ocx().Tags._resource._schema['properties']['name']['maxLength'] + 1
    invalid_tag = 'a' * too_many_characters
    initial_appendables = { 'metadata': {}, 'tags': [invalid_tag], 'valid_tags': [], 'custom_metadata': {}, 'valid_metadata': {}}
    with pytest.raises(ValidationError):
        validate_tags(initial_appendables, ocx())


def test_validate_metadata_valid():
    initial_appendables = { 'metadata': {'starred': 'true'}, 'tags': [], 'valid_tags': [], 'custom_metadata': {}, 'valid_metadata': {}}
    final_appendables = dict(initial_appendables)
    final_appendables['custom_metadata']['starred'] = True
    validate_metadata(initial_appendables, ocx())
    assert initial_appendables == final_appendables


def test_validate_metadata_not_present():
    initial_appendables = { 'tags': [], 'valid_tags': [], 'custom_metadata': {}, 'valid_metadata': {}}
    final_appendables = dict(initial_appendables)
    validate_metadata(initial_appendables, ocx())
    assert initial_appendables == final_appendables


def test_validate_metadata_blacklisted():
    initial_appendables = { 'metadata': {'$uri': 'invalid_entry'}, 'tags': [], 'valid_tags': [], 'custom_metadata': {}, 'valid_metadata': {}}
    with pytest.raises(ValidationError):
        validate_metadata(initial_appendables, ocx())


def test_validate_metadata_custom_valid():
    initial_appendables = { 'metadata': {'foo': 'bar'}, 'tags': [], 'valid_tags': [], 'custom_metadata': {}, 'valid_metadata': {}}
    final_appendables = dict(initial_appendables)
    final_appendables['custom_metadata']['foo'] ='bar'
    validate_metadata(initial_appendables, ocx())
    assert initial_appendables == final_appendables


def test_validate_metadata_against_schema():
    schema_props = ocx().Metadata._resource._schema['properties']
    with mock.patch('onecodex.metadata_upload.validate_enum') as mock_enum:
        validate_metadata_against_schema(schema_props, 'platform', 'Illumina')
        assert mock_enum.call_count == 1

    with mock.patch('onecodex.metadata_upload.validate_number') as mock_number:
        validate_metadata_against_schema(schema_props, 'location_lat', '90')
        assert mock_number.call_count == 1

    with mock.patch('onecodex.metadata_upload.validate_boolean') as mock_bool:
        validate_metadata_against_schema(schema_props, 'starred', 'true')
        assert mock_bool.call_count == 1

    with mock.patch('onecodex.metadata_upload.validate_datetime') as mock_datetime:
        validate_metadata_against_schema(schema_props, 'date_collected', '2018,05,08,12,12')
        assert mock_datetime.call_count == 1
    string_response = validate_metadata_against_schema(schema_props, 'type', 'foo')
    assert string_response == 'foo'

