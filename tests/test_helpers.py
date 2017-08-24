import pandas as pd
import pytest

from onecodex.exceptions import OneCodexException
from onecodex.helpers import collate_classification_results, normalize_classifications
from onecodex.models import Classifications


def test_normalize_classifications(ocx, api_data):
    classifications = [ocx.Classifications.get('45a573fb7833449a'),
                       ocx.Classifications.get('593601a797914cbf')]
    normed_classifications, metadata = normalize_classifications(classifications)
    assert len(normed_classifications) == 2
    assert len(metadata) == 2
    assert isinstance(metadata, pd.DataFrame)
    assert (metadata['_display_name'] == metadata['name']).all()

    for n in normed_classifications:
        assert isinstance(n, Classifications)


def test_normalize_analyses_and_samples(ocx, api_data):
    analyses = [ocx.Analyses.get('45a573fb7833449a'),
                ocx.Analyses.get('593601a797914cbf')]
    normed_classifications, metadata = normalize_classifications(analyses)
    for n in normed_classifications:
        assert isinstance(n, Classifications)

    normed_classifications, metadata = normalize_classifications([a.sample for a in analyses])
    for n in normed_classifications:
        assert isinstance(n, Classifications)


def test_normalize_classifications_labeling(ocx, api_data):
    cs = [ocx.Classifications.get('45a573fb7833449a'),
          ocx.Classifications.get('593601a797914cbf')]
    _, metadata = normalize_classifications(cs, label='created_at')
    assert (metadata['_display_name'] == metadata['created_at']).all()

    with pytest.raises(OneCodexException):
        _, metadata = normalize_classifications(cs, label='nonexisting')

    _, metadata = normalize_classifications(cs, label=lambda x: x.sample.filename)
    assert (metadata['_display_name'] == [c.sample.filename for c in cs]).all()


def test_unsuccessful_classification(ocx, api_data):
    classifications = [ocx.Classifications.get('45a573fb7833449a'),
                       ocx.Classifications.get('593601a797914cbf')]
    classifications[0]._resource.success = False
    with pytest.warns(UserWarning):
        normed_classifications, _ = normalize_classifications(classifications)
        assert len(normed_classifications) == 1

    normed_classifications, _ = normalize_classifications(classifications, skip_missing=False)
    assert len(normed_classifications) == 2


@pytest.mark.parametrize('rank,field,remove_zeros,val1,val2,val3,val4', [
    (None, 'readcount_w_children', True, [3], [80], [3], [80]),
    (None, 'readcount', True, [], [], [3], [80]),
    ('genus', 'readcount_w_children', True, [3], [80], [], []),
    ('genus', 'readcount', True, [], [], [], []),
    ('genus', 'readcount', False, [0], [0], [], []),
])
def test_collate_classifications(ocx, api_data, rank, field, remove_zeros,
                                 val1, val2, val3, val4):
    classifications = [ocx.Classifications.get('45a573fb7833449a'),
                       ocx.Classifications.get('593601a797914cbf')]
    df, tax_info = collate_classification_results(classifications, rank=rank, field=field,
                                                  remove_zeros=remove_zeros)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 2
    assert df.columns.names == ['tax_id']  # Don't bother w/ a multi-index by default

    SAMPLE_IDS = ['45a573fb7833449a', '593601a797914cbf', '45a573fb7833449a', '593601a797914cbf']
    TAX_IDS = ['1279', '1279', '1078083', '1078083']
    for ix, val in enumerate([val1, val2, val3, val4]):
        if val:
            assert df.loc[SAMPLE_IDS[ix], TAX_IDS[ix]] == val
        else:
            with pytest.raises(KeyError):
                df.loc[SAMPLE_IDS[ix], TAX_IDS[ix]]
