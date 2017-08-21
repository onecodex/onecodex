import pytest
import skbio

from onecodex.distance import (alpha_counts, beta_counts, beta_diversity,
                               chao1, simpson, unifrac, jaccard, braycurtis, cityblock)


def test_alpha_counts(ocx, api_data):
    analysis = ocx.Classifications.get('45a573fb7833449a')
    counts, ids = alpha_counts(analysis)
    assert counts == [3]
    assert ids == [analysis.sample.metadata.name]

    counts, ids = alpha_counts(analysis, rank='genus')
    assert counts == [3]
    counts, ids = alpha_counts(analysis, rank='genus', field='readcount')
    assert counts == [0]

    counts, ids = alpha_counts(ocx.Classifications.get('593601a797914cbf'))
    assert counts == [80]


def test_beta_counts(ocx, api_data):
    analyses = [ocx.Classifications.get('45a573fb7833449a'),
                ocx.Classifications.get('593601a797914cbf')]

    vectors, tax_ids, ids = beta_counts(analyses)
    assert vectors == [[3], [80]]
    assert tax_ids == ['1078083']

    vectors, tax_ids, ids = beta_counts(analyses, rank='genus')
    assert vectors == [[3], [80]]
    assert tax_ids == ['1279']

    vectors, tax_ids, ids = beta_counts(analyses, rank='genus', field='readcount')
    assert vectors == [[], []]
    assert tax_ids == []
    assert list(ids) == [analyses[0].id, analyses[1].id]


def test_beta_diversity(ocx, api_data):
    analyses = [ocx.Classifications.get('45a573fb7833449a'),
                ocx.Classifications.get('593601a797914cbf')]
    dm = beta_diversity(analyses, 'braycurtis')
    assert isinstance(dm, skbio.stats.distance._base.DistanceMatrix)
    assert list(dm.data[0, :]) == [0.0, 0.92771084337349397]

    with pytest.raises(ValueError):
        beta_diversity(analyses, 'invalid_metric')


@pytest.mark.parametrize('f,value', [
    (chao1, 1.0),
    (simpson, 0.0)
])
def test_alphas(ocx, api_data, f, value):
    analysis = ocx.Classifications.get('45a573fb7833449a')
    assert f(analysis) == value


@pytest.mark.parametrize('f,value,kwargs', [
    (braycurtis, [0.0, 0.92771084337349397], {}),
    (cityblock, [0.0, 77.0], {}),
    (jaccard, [0.0, 1.0], {}),
    (unifrac, [0.0, 0.0], {'strict': False}),
])
def test_unifrac(ocx, api_data, f, value, kwargs):
    analyses = [ocx.Classifications.get('45a573fb7833449a'),
                ocx.Classifications.get('593601a797914cbf')]
    dm = f(analyses, **kwargs)
    assert isinstance(dm, skbio.stats.distance._base.DistanceMatrix)
    assert list(dm._data[0, :]) == value
