import pytest
import warnings

from onecodex.exceptions import OneCodexException
from onecodex.viz import plot_heatmap, plot_metadata, plot_pca


# TODO: somehow check the graphical output (JSON output to Vega?)
def test_plot_metadata(ocx, api_data):
    analyses = [ocx.Classifications.get('6579e99943f84ad2'),
                ocx.Classifications.get('b50c176668234fe7'),
                ocx.Classifications.get('e0422602de41479f')]

    samples = [a.sample for a in analyses]

    # expect a warning (issue only once) when specifying alphadiv without passing rank
    with warnings.catch_warnings(record=True) as w:
        plot_metadata(analyses, alphadiv='simpson')
        assert len(w) == 1
        assert 'using species' in str(w[-1].message)

    # should resolve samples to classifications via normalize_classifications()
    plot_metadata(samples, alphadiv='simpson')

    # try time, boolean, and numerical types for x-axis
    plot_metadata(samples, x='date_sequenced', alphadiv='chao1')
    plot_metadata(samples, x='starred', alphadiv='chao1')
    plot_metadata(samples, x='totalige', alphadiv='chao1')

    # taxid and taxon on vertical axis
    plot_metadata(samples, taxid=1279)  # should coerce to string internally
    plot_metadata(samples, taxid='1279')
    plot_metadata(samples, taxon='Staphylococcus')

    # force a scatter to be a boxplot and vice versa
    plot_metadata(samples, x='totalige', alphadiv='chao1', boxplot=True)
    plot_metadata(samples, x='starred', alphadiv='chao1', scatter=True)

    # test plot labeling
    plot_metadata(
        samples, alphadiv='simpson', title='my title', xlabel='my xlabel', ylabel='my ylabel'
    )


def test_plot_metadata_exceptions(ocx, api_data):
    analyses = [ocx.Classifications.get('6579e99943f84ad2'),
                ocx.Classifications.get('b50c176668234fe7'),
                ocx.Classifications.get('e0422602de41479f')]

    # cant specify more than one vert axis parameter simultaneously
    with pytest.raises(OneCodexException):
        plot_metadata(analyses, alphadiv='simpson', metadata='_display_name')

    with pytest.raises(OneCodexException):
        plot_metadata(analyses, taxid='1279', taxon='Staphylococcus')

    # cant force boxplot and scatter plot simultaneously
    with pytest.raises(OneCodexException):
        plot_metadata(analyses, taxid='1279', boxplot=True, scatter=True)

    # alphadiv metric must be simpson or chao1
    with pytest.raises(OneCodexException):
        plot_metadata(analyses, alphadiv='does_not_exist')

    # vert axis metadata must exist and be numeric
    with pytest.raises(OneCodexException) as e:
        plot_metadata(analyses, metadata='does_not_exist')
    assert 'not in metadata' in str(e.value)

    with pytest.raises(OneCodexException) as e:
        plot_metadata(analyses, metadata='delivery')
    assert 'must be numerical' in str(e.value)

    # taxid and taxon that don't exist
    with pytest.raises(OneCodexException) as e:
        plot_metadata(analyses, taxid='does_not_exist')
    assert 'not found in analyses' in str(e.value)

    with pytest.raises(OneCodexException) as e:
        plot_metadata(analyses, taxon='does_not_exist')
    assert 'not found in analyses' in str(e.value)

    # taxid and taxon that don't exist in the provided rank
    with pytest.raises(OneCodexException) as e:
        plot_metadata(analyses, taxid='1279', rank='family')
    assert 'not found in analyses' in str(e.value)

    with pytest.raises(OneCodexException) as e:
        plot_metadata(analyses, taxon='Staphylococcus', rank='family')
    assert 'not found in analyses' in str(e.value)

    # single x-axis param must exist
    with pytest.raises(OneCodexException) as e:
        plot_metadata(analyses, x='does_not_exist', alphadiv='simpson')
    assert 'not in metadata' in str(e.value)

    # when multiple x-axis params, all must exist and be categorical
    with pytest.raises(OneCodexException) as e:
        plot_metadata(analyses, x=['vegetables', 'does_not_exist'], alphadiv='simpson')
    assert 'not in metadata' in str(e.value)

    with pytest.raises(OneCodexException) as e:
        plot_metadata(analyses, x=['eggs', 'totalige'], alphadiv='simpson')
    assert 'must be categorical' in str(e.value)


def test_plot_pca(ocx, api_data):
    # plot_pca requires >1 classification result and >1 taxon within the specified rank
    analyses = [ocx.Classifications.get('6579e99943f84ad2'),
                ocx.Classifications.get('b50c176668234fe7'),
                ocx.Classifications.get('e0422602de41479f')]

    plot_pca(analyses)

    # test plot labeling
    plot_pca(analyses, title='my title', xlabel='my xlabel', ylabel='my ylabel')

    # test changing size/colors by metadata
    plot_pca(analyses, color='geo_loc_name', size='totalige')

    # test changing size/colors by taxid
    plot_pca(analyses, color='taxid_1279', size='taxid_816')


def test_plot_pca_exceptions(ocx, api_data):
    # large, full metadata analyses
    analyses = [ocx.Classifications.get('6579e99943f84ad2'),
                ocx.Classifications.get('b50c176668234fe7'),
                ocx.Classifications.get('e0422602de41479f')]

    # tiny analyses that aren't plottable
    tiny_analyses = [ocx.Classifications.get('45a573fb7833449a'),
                     ocx.Classifications.get('593601a797914cbf')]

    # taxonomic rank must be specified
    with pytest.raises(OneCodexException) as e:
        plot_pca(analyses, rank=None)
    assert 'specify a taxonomic rank' in str(e.value)

    # must have at least two samples
    with pytest.raises(OneCodexException) as e:
        plot_pca(analyses[0])
    assert 'requires 2 or more' in str(e.value)

    with pytest.raises(OneCodexException) as e:
        plot_pca([analyses[0]])
    assert 'requires 2 or more' in str(e.value)

    # samples must have at least two taxa at this rank
    with pytest.raises(OneCodexException) as e:
        plot_pca(tiny_analyses, rank='species')
    assert 'at least 2 for PCA' in str(e.value)

    # color/size/tooltips with invalid metadata fields or taxids
    for k in ('color', 'size', 'tooltip'):
        kwargs = {k: 'does_not_exist'}
        with pytest.raises(OneCodexException) as e:
            plot_pca(analyses, **kwargs)
        assert 'not found in metadata' in str(e.value)

        kwargs = {k: 'taxid_does_not_exist'}
        with pytest.raises(OneCodexException) as e:
            plot_pca(analyses, **kwargs)
        assert 'not found in analyses' in str(e.value)

    # specifying taxid at incorrect rank relative to what's been specified
    with pytest.raises(OneCodexException) as e:
        plot_pca(analyses, rank='family', color='taxid_1279', size='taxid_816')
    assert 'not found in analyses' in str(e.value)


def test_plot_heatmap(ocx, api_data):
    analyses = [ocx.Classifications.get('6579e99943f84ad2'),
                ocx.Classifications.get('b50c176668234fe7'),
                ocx.Classifications.get('e0422602de41479f')]

    plot_heatmap(analyses, top_n=10, threshold=None)
    plot_heatmap(analyses, top_n=None, threshold=0.1)
    plot_heatmap(analyses, top_n=10, threshold=0.1)
    plot_heatmap(analyses, title='my title', xlabel='my xlabel', ylabel='my ylabel')


def test_plot_heatmap_exceptions(ocx, api_data):
    analyses = [ocx.Classifications.get('6579e99943f84ad2'),
                ocx.Classifications.get('b50c176668234fe7'),
                ocx.Classifications.get('e0422602de41479f')]

    # taxonomic rank must be specified
    with pytest.raises(OneCodexException) as e:
        plot_heatmap(analyses, rank=None)
    assert 'specify a taxonomic rank' in str(e.value)

    # must have at least two samples
    with pytest.raises(OneCodexException) as e:
        plot_heatmap(analyses[0])
    assert 'requires 2 or more' in str(e.value)

    with pytest.raises(OneCodexException) as e:
        plot_heatmap([analyses[0]])
    assert 'requires 2 or more' in str(e.value)

    # must specify at least threshold or top_n
    with pytest.raises(OneCodexException) as e:
        plot_heatmap(analyses, top_n=None, threshold=None)
    assert 'specify one of' in str(e.value)

    # tooltip with invalid metadata fields or taxids
    with pytest.raises(OneCodexException) as e:
        plot_heatmap(analyses, tooltip='does_not_exist')
    assert 'not found in metadata' in str(e.value)

    with pytest.raises(NotImplementedError) as e:
        plot_heatmap(analyses, tooltip='taxid_does_not_exist')
    assert 'not supported' in str(e.value)
