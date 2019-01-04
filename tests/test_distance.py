import pandas as pd
import pytest
import skbio

from onecodex.exceptions import OneCodexException


@pytest.mark.parametrize('metric,value,kwargs', [
    ('simpson', [0.7887047744107811, 0.8459005545543041, 0.3722838525151998], {}),
    ('chao1', [402.5833333333333, 374.2830188679245, 308.10714285714283], {}),
])
def test_alpha_diversity(ocx_w_enhanced, api_data, metric, value, kwargs):
    samples = ocx_w_enhanced.Samples.where(project='4b53797444f846c4')
    divs = samples.alpha_diversity(metric=metric, **kwargs)
    assert isinstance(divs, pd.DataFrame)
    assert divs[metric].tolist() == value


def test_alpha_diversity_exceptions(ocx_w_enhanced, api_data):
    samples = ocx_w_enhanced.Samples.where(project='4b53797444f846c4')

    # should fail if data has been normalized
    with pytest.raises(OneCodexException) as e:
        samples.results(normalize=True).ocx.alpha_diversity('simpson')
    assert 'requires unnormalized' in str(e.value)

    # must be a metric that exists
    with pytest.raises(OneCodexException) as e:
        samples.alpha_diversity('does_not_exist')
    assert 'metric must be one of' in str(e.value)


@pytest.mark.parametrize('metric,value,kwargs', [
    ('braycurtis', [0.758579937018798, 0.46261493509445945, 0.7603369765359447], {}),
    ('cityblock', [18680367.0, 13918888.0, 16237777.0], {}),
    ('jaccard', [0.9823788546255506, 0.9638242894056848, 0.9678217821782178], {}),
])
def test_beta_diversity(ocx_w_enhanced, api_data, metric, value, kwargs):
    samples = ocx_w_enhanced.Samples.where(project='4b53797444f846c4')
    dm = samples.beta_diversity(metric=metric, **kwargs)
    assert isinstance(dm, skbio.stats.distance._base.DistanceMatrix)
    assert dm.condensed_form().tolist() == value


def test_beta_diversity_exceptions(ocx_w_enhanced, api_data):
    samples = ocx_w_enhanced.Samples.where(project='4b53797444f846c4')

    # should fail if data has been normalized
    with pytest.raises(OneCodexException) as e:
        samples.results(normalize=True).ocx.beta_diversity('braycurtis')
    assert 'requires unnormalized' in str(e.value)

    # must be a metric that exists
    with pytest.raises(OneCodexException) as e:
        samples.beta_diversity('does_not_exist')
    assert 'metric must be one of' in str(e.value)


@pytest.mark.parametrize('weighted,value', [
    (False, [0.380050505050505, 0.3538011695906433, 0.4122448979591837]),
    (True, [6.418015325847228, 5.209385639879688, 8.104451028880142]),
])
def test_unifrac(ocx_w_enhanced, api_data, value, weighted):
    samples = ocx_w_enhanced.Samples.where(project='4b53797444f846c4')
    dm = samples.unifrac(weighted=weighted)
    assert isinstance(dm, skbio.stats.distance._base.DistanceMatrix)
    assert dm.condensed_form().tolist() == value


def test_unifrac_exceptions(ocx_w_enhanced, api_data):
    samples = ocx_w_enhanced.Samples.where(project='4b53797444f846c4')

    # should fail if data has been normalized
    with pytest.raises(OneCodexException) as e:
        samples.results(normalize=True).ocx.unifrac()
    assert 'requires unnormalized' in str(e.value)
