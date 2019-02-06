from click import BadParameter
from functools import partial
import pytest

from onecodex.utils import snake_case, check_for_allowed_file, valid_api_key


def test_check_allowed_file():
    # bad ones
    with pytest.raises(SystemExit):
        check_for_allowed_file("file.bam")
        check_for_allowed_file("file")

    # good ones
    check_for_allowed_file("file.fastq")
    check_for_allowed_file("file.fastq.gz")


def test_is_valid_api_key():
    empty_key = ""
    short_key = "123"
    long_key = "123abc123abc123abc123abc123abc123abc123abc123abc123abc"
    good_key = "123abc123abc123abc123abc123abc32"

    # its a click callback so it expects some other params
    valid_api_key_partial = partial(valid_api_key, None, None)

    for key in [empty_key, short_key, long_key]:
        with pytest.raises(BadParameter):
            valid_api_key_partial(key)

    assert good_key == valid_api_key_partial(good_key)


@pytest.mark.parametrize('resource,uris', [
    ('Samples', []),
    ('Samples', ['761bc54b97f64980']),
    ('Analyses', []),
    ('Analyses', ['45a573fb7833449a']),
    ('Markerpanels', [])
])
def test_fetcher(ocx, api_data, resource, uris):
    if len(uris) == 0:
        pass
    else:
        for uri in uris:
            resource_class = getattr(ocx, resource)
            instance = resource_class.get(uri)
            assert instance is not None


def test_snake_case():
    test_cases = ['SnakeCase', 'snakeCase', 'SNAKE_CASE']
    for test_case in test_cases:
        assert snake_case(test_case) == 'snake_case'
