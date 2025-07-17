import pytest
import responses


@pytest.fixture
def sample_response_data():
    """Sample response data that would come from ?expand=all."""
    return {
        "GET::api/v1/samples/0ee172af60e84f61": {
            "$uri": "/api/v1/samples/0ee172af60e84f61",
            "created_at": "2023-01-01T00:00:00+00:00",
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
            "status": "available",
        }
    }


@pytest.fixture
def detailed_sample_data():
    """More detailed sample data for testing."""
    return {
        "GET::api/v1/samples/0ee172af60e84f61": {
            "$uri": "/api/v1/samples/0ee172af60e84f61",
            "created_at": "2023-01-01T00:00:00+00:00",
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
            "status": "available",
        }
    }


def test_sample_metadata_sample_no_additional_requests(
    ocx, custom_mock_requests, sample_response_data
):
    """Test that getting a sample with expand=all doesn't trigger additional requests."""
    with custom_mock_requests(sample_response_data):
        sample = ocx.Samples.get("0ee172af60e84f61")

        # Access metadata - this should be already resolved from expand=all
        metadata = sample.metadata
        assert metadata is not None
        assert metadata.platform == "Illumina"
        assert metadata.field_uri == "/api/v1/metadata/dd0aa4fdd11541d1"

        # Access the circular reference - should resolve to the original sample
        # without making additional requests
        circular_sample = sample.metadata.sample

        # The circular sample should be the same as the original sample,
        # and populated from a single expanded request
        assert circular_sample is sample
        assert circular_sample.field_uri == sample.field_uri
        assert len(responses.calls) == 1
        assert responses.calls[0].request.url is not None
        assert "?expand=all" in responses.calls[0].request.url


def test_sample_metadata_sample_with_detailed_data(ocx, custom_mock_requests, detailed_sample_data):
    """Test with more detailed mock data structure."""
    with custom_mock_requests(detailed_sample_data):
        # Get sample and verify behavior
        sample = ocx.Samples.get("0ee172af60e84f61")

        # Should have made exactly one request
        assert len(responses.calls) == 1

        # Access metadata - should be resolved object, not ApiRef
        metadata = sample.metadata
        assert hasattr(metadata, "field_uri")
        assert metadata.field_uri == "/api/v1/metadata/dd0aa4fdd11541d1"

        # Still only one request
        assert len(responses.calls) == 1

        # Access circular reference
        circular_sample = sample.metadata.sample

        # Should still be only one request (no additional API call for circular ref)
        assert len(responses.calls) == 1

        # Should be the same object
        assert circular_sample is sample


def test_apiref_circular_resolution_detection(ocx, custom_mock_requests):
    """Test that circular references are detected during lazy loading."""
    # Mock response for sample that has lazy-loaded metadata
    unexpanded_response_data = {
        "GET::api/v1/samples/0ee172af60e84f61": {
            "$uri": "/api/v1/samples/0ee172af60e84f61",
            "created_at": "2023-01-01T00:00:00+00:00",
            "filename": "test.fastq",
            "metadata": {"$ref": "/api/v1/metadata/dd0aa4fdd11541d1"},
            "owner": {"$ref": "/api/v1/users/5891ee65711c4d5e"},
            "primary_classification": {"$ref": "/api/v1/classifications/cb18d91b82214e11"},
            "project": None,
            "tags": [],
            "visibility": "private",
            "size": 1000,
            "status": "available",
        },
        "GET::api/v1/metadata/dd0aa4fdd11541d1": {
            "$uri": "/api/v1/metadata/dd0aa4fdd11541d1",
            "sample": {"$ref": "/api/v1/samples/0ee172af60e84f61"},
            "platform": "Illumina",
        },
    }

    with custom_mock_requests(unexpanded_response_data):
        # Get the sample - this triggers lazy loading
        sample = ocx.Samples.get("0ee172af60e84f61")
        assert len(responses.calls) == 1

        # Access metadata - should trigger one more API call
        metadata = sample.metadata
        assert len(responses.calls) == 2

        # Access the circular reference - should NOT trigger another API call
        # because it should detect the circular reference and return the original sample
        circular_sample = metadata.sample

        # Should still be only two requests (no additional call for circular ref)
        assert len(responses.calls) == 2

        # Should be the same object
        assert circular_sample is sample


def test_multiple_circular_accesses_use_cache(ocx, custom_mock_requests, sample_response_data):
    """Test that multiple accesses to the same circular reference use caching."""

    with custom_mock_requests(sample_response_data):
        sample = ocx.Samples.get("0ee172af60e84f61")
        assert len(responses.calls) == 1

        # Access the circular reference multiple times
        ref1 = sample.metadata.sample
        ref2 = sample.metadata.sample
        ref3 = sample.metadata.sample

        # Should still be only one request
        assert len(responses.calls) == 1

        # All references should be the same object
        assert ref1 is sample
        assert ref2 is sample
        assert ref3 is sample
        assert ref1 is ref2 is ref3
