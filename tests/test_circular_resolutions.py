import pytest
from unittest.mock import MagicMock
from onecodex import Api


@pytest.fixture
def mock_api():
    """Create a mock API object."""
    mock_api = MagicMock()
    mock_api._base_url = "https://app.onecodex.com"
    return mock_api


@pytest.fixture
def mock_client():
    """Create a mock HTTP client."""
    return MagicMock()


@pytest.fixture
def sample_response_data():
    """Sample response data that would come from ?expand=all."""
    return {
        "$uri": "/api/v1/samples/0ee172af60e84f61",
        "created_at": "2023-01-01T00:00:00Z",
        "filename": "test.fastq",
        "metadata": {
            "$uri": "/api/v1/metadata/dd0aa4fdd11541d1",
            "sample": {"$ref": "/api/v1/samples/0ee172af60e84f61"},  # Circular ref as $ref
            "platform": "Illumina",
        },
        "owner": {"$ref": "/api/v1/users/5891ee65711c4d5e"},
        "primary_classification": {"$ref": "/api/v1/classifications/cb18d91b82214e11"},
        "project": None,
        "tags": [],
        "visibility": "private",
        "size": 1000,
    }


@pytest.fixture
def detailed_sample_data():
    """More detailed sample data for testing."""
    return {
        "$uri": "/api/v1/samples/0ee172af60e84f61",
        "created_at": "2023-01-01T00:00:00Z",
        "filename": "test.fastq",
        "metadata": {
            "$uri": "/api/v1/metadata/dd0aa4fdd11541d1",
            "sample": {"$ref": "/api/v1/samples/0ee172af60e84f61"},
            "platform": "Illumina",
            "library_type": "WGS",
            "sample_type": "Metagenomic",  # Use a valid enum value
        },
        "owner": {"$ref": "/api/v1/users/5891ee65711c4d5e"},
        "primary_classification": {"$ref": "/api/v1/classifications/cb18d91b82214e11"},
        "project": None,
        "tags": [],
        "visibility": "private",
        "size": 1000,
    }


@pytest.fixture
def configured_api(mock_api, mock_client):
    """Create a configured API instance with mocked components."""
    ocx = Api()
    ocx.Samples._client = mock_client
    ocx.Samples._api = mock_api
    ocx.Metadata._client = mock_client
    ocx.Metadata._api = mock_api
    return ocx


def test_sample_metadata_sample_no_additional_requests(
    configured_api, mock_client, sample_response_data
):
    """Test that sample.metadata.sample doesn't trigger additional API requests."""

    # Setup the mock response
    mock_response = MagicMock()
    mock_response.json.return_value = sample_response_data
    mock_client.get.return_value = mock_response

    # Get the sample
    sample = configured_api.Samples.get("0ee172af60e84f61")

    # Verify the initial GET request was made
    assert mock_client.get.call_count == 1
    mock_client.get.assert_called_with(
        "https://app.onecodex.com/api/v1/samples/0ee172af60e84f61?expand=all"
    )

    # Access metadata - this should be already resolved from expand=all
    metadata = sample.metadata
    assert metadata is not None
    assert hasattr(metadata, "field_uri")
    assert metadata.field_uri == "/api/v1/metadata/dd0aa4fdd11541d1"

    # No additional requests should have been made
    assert mock_client.get.call_count == 1

    # Access the circular reference - this should resolve to the original sample
    # without making additional requests
    assert id(sample) == id(sample.metadata.sample)
    circular_sample = sample.metadata.sample

    # Verify no additional requests were made
    assert mock_client.get.call_count == 1

    # The circular sample should be the same as the original sample
    assert circular_sample is sample
    assert circular_sample.field_uri == sample.field_uri


def test_sample_metadata_sample_with_detailed_data(
    configured_api, mock_client, detailed_sample_data
):
    """Test with more detailed mock data structure."""

    # Setup mock response
    mock_response = MagicMock()
    mock_response.json.return_value = detailed_sample_data
    mock_client.get.return_value = mock_response

    # Get sample and verify behavior
    sample = configured_api.Samples.get("0ee172af60e84f61")

    # Should have made exactly one request
    assert mock_client.get.call_count == 1

    # Access metadata - should be resolved object, not ApiRef
    metadata = sample.metadata
    assert hasattr(metadata, "field_uri")
    assert metadata.field_uri == "/api/v1/metadata/dd0aa4fdd11541d1"

    # Still only one request
    assert mock_client.get.call_count == 1

    # Access circular reference
    circular_sample = sample.metadata.sample

    # Should still be only one request (no additional API call for circular ref)
    assert mock_client.get.call_count == 1

    # Should be the same object
    assert circular_sample is sample


def test_apiref_circular_resolution_detection(configured_api, mock_client):
    """Test that circular references are detected during lazy loading."""
    # Mock response for metadata that has a lazy-loaded sample reference
    metadata_response = {
        "$uri": "/api/v1/metadata/dd0aa4fdd11541d1",
        "sample": {"$ref": "/api/v1/samples/0ee172af60e84f61"},
        "platform": "Illumina",
    }

    # Mock response for sample that has lazy-loaded metadata
    sample_response = {
        "$uri": "/api/v1/samples/0ee172af60e84f61",
        "created_at": "2023-01-01T00:00:00Z",
        "filename": "test.fastq",
        "metadata": {"$ref": "/api/v1/metadata/dd0aa4fdd11541d1"},
        "owner": {"$ref": "/api/v1/users/5891ee65711c4d5e"},
        "primary_classification": {"$ref": "/api/v1/classifications/cb18d91b82214e11"},
        "project": None,
        "tags": [],
        "visibility": "private",
        "size": 1000,
    }

    def mock_get_response(url):
        mock_response = MagicMock()
        if "metadata" in url:
            mock_response.json.return_value = metadata_response
        else:
            mock_response.json.return_value = sample_response
        return mock_response

    mock_client.get.side_effect = mock_get_response

    # Get the sample - this triggers lazy loading
    sample = configured_api.Samples.get("0ee172af60e84f61")
    assert mock_client.get.call_count == 1

    # Access metadata - should trigger one more API call
    metadata = sample.metadata
    assert mock_client.get.call_count == 2

    # Access the circular reference - should NOT trigger another API call
    # because it should detect the circular reference and return the original sample
    circular_sample = metadata.sample

    # Should still be only two requests (no additional call for circular ref)
    assert mock_client.get.call_count == 2

    # Should be the same object
    assert circular_sample is sample


def test_multiple_circular_accesses_use_cache(configured_api, mock_client, sample_response_data):
    """Test that multiple accesses to the same circular reference use caching."""

    # Setup the mock response
    mock_response = MagicMock()
    mock_response.json.return_value = sample_response_data
    mock_client.get.return_value = mock_response

    # Get the sample
    sample = configured_api.Samples.get("0ee172af60e84f61")
    assert mock_client.get.call_count == 1

    # Access the circular reference multiple times
    ref1 = sample.metadata.sample
    ref2 = sample.metadata.sample
    ref3 = sample.metadata.sample

    # Should still be only one request
    assert mock_client.get.call_count == 1

    # All references should be the same object
    assert ref1 is sample
    assert ref2 is sample
    assert ref3 is sample
    assert ref1 is ref2 is ref3
