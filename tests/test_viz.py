import pytest; pytest.importorskip('pandas')  # noqa

from onecodex.exceptions import OneCodexException


# TODO: somehow check the graphical output (JSON output to Vega?)
def test_plot_metadata(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')

    # should resolve samples to classifications via normalize_classifications()
    samples.plot_metadata(vaxis='simpson')

    # try time, boolean, and numerical types for x-axis
    samples.plot_metadata(haxis='date_sequenced', vaxis='chao1')
    samples.plot_metadata(haxis='starred', vaxis='chao1')
    samples.plot_metadata(haxis='totalige', vaxis='chao1')

    # taxid and taxon on vertical axis
    samples.plot_metadata(vaxis=1279)  # should coerce to string internally
    samples.plot_metadata(vaxis='1279')
    samples.plot_metadata(vaxis='Staphylococcus')

    # force a scatter to be a boxplot and vice versa
    samples.plot_metadata(haxis='totalige', vaxis='chao1', plot_type='boxplot')
    samples.plot_metadata(haxis='starred', vaxis='chao1', plot_type='scatter')

    # test plot labeling
    samples.plot_metadata(
        vaxis='simpson', title='my title', xlabel='my xlabel', ylabel='my ylabel'
    )


def test_plot_metadata_exceptions(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')

    # expect error if rank is None, since that could lead to weird results
    with pytest.raises(OneCodexException) as e:
        samples.plot_metadata(vaxis='simpson', rank=None)
    assert 'specify a rank' in str(e.value)

    # vert axis does not exist
    with pytest.raises(OneCodexException):
        samples.plot_metadata(vaxis='does_not_exist')

    # horiz axis does not exist
    with pytest.raises(OneCodexException):
        samples.plot_metadata(haxis='does_not_exist')


def test_plot_pca(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')

    samples.plot_pca()

    # test plot labeling
    samples.plot_pca(title='my title', xlabel='my xlabel', ylabel='my ylabel')

    # test changing size/colors by metadata
    samples.plot_pca(color='geo_loc_name', size='totalige')

    # test changing size/colors by tax id
    samples.plot_pca(color='1279', size='816')
    samples.plot_pca(color=1279, size=816)

    # test changing size/colors by taxon name
    samples.plot_pca(color='Bacteroides', size='Prevotella')

    # test tooltips
    samples.plot_pca(tooltip=['totalige', 'vegetables', 'Prevotella', '816'])


def test_plot_pca_exceptions(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')

    # expect error if rank is None, since that could lead to weird results
    with pytest.raises(OneCodexException) as e:
        samples.plot_pca(rank=None)
    assert 'specify a rank' in str(e.value)

    # must have at least two samples
    with pytest.raises(OneCodexException) as e:
        samples[:1].plot_pca()
    assert 'requires 2 or more' in str(e.value)

    # samples must have at least two taxa at this rank
    with pytest.raises(OneCodexException) as e:
        samples.to_df(top_n=1).ocx.plot_pca()
    assert 'at least 2 for PCA' in str(e.value)

    # color/size/tooltips with invalid metadata fields or taxids
    for k in ('color', 'size', 'tooltip'):
        kwargs = {k: 'does_not_exist'}
        with pytest.raises(OneCodexException) as e:
            samples.plot_pca(**kwargs)
        assert 'not found' in str(e.value)

        kwargs = {k: '487527863'}
        with pytest.raises(OneCodexException) as e:
            samples.plot_pca(**kwargs)
        assert 'not found' in str(e.value)


def test_plot_heatmap(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')

    samples.plot_heatmap(top_n=10, threshold=None)
    samples.plot_heatmap(top_n=None, threshold=0.1)
    samples.plot_heatmap(top_n=10, threshold=0.1)
    samples.plot_heatmap(top_n=10, haxis='eggs')
    samples.plot_heatmap(title='my title', xlabel='my xlabel', ylabel='my ylabel')


def test_plot_heatmap_exceptions(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')

    # expect error if rank is None, since that could lead to weird results
    with pytest.raises(OneCodexException) as e:
        samples.plot_heatmap(rank=None)
    assert 'specify a rank' in str(e.value)

    # must have at least two samples
    with pytest.raises(OneCodexException) as e:
        samples[:1].plot_heatmap()
    assert 'requires 2 or more' in str(e.value)

    # must specify at least threshold or top_n
    with pytest.raises(OneCodexException) as e:
        samples.plot_heatmap(top_n=None, threshold=None)
    assert 'specify at least one of' in str(e.value)

    # tooltip with invalid metadata fields or taxids
    with pytest.raises(OneCodexException) as e:
        samples.plot_heatmap(tooltip='does_not_exist')
    assert 'not found' in str(e.value)


@pytest.mark.parametrize('metric', ('braycurtis', 'jaccard', 'unifrac', 'manhattan'))
def test_plot_distance(ocx, api_data, metric):
    samples = ocx.Samples.where(project='4b53797444f846c4')

    samples.plot_distance(metric=metric, xlabel='my xlabel', ylabel='my ylabel', title='my title')


def test_plot_distance_exceptions(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')

    # expect error if rank is None, since that could lead to weird results
    with pytest.raises(OneCodexException) as e:
        samples.plot_distance(rank=None)
    assert 'specify a rank' in str(e.value)

    # only some metrics allowed
    with pytest.raises(OneCodexException) as e:
        samples.plot_distance(metric='simpson')
    assert 'must be one of' in str(e.value)

    # need more than one analysis
    with pytest.raises(OneCodexException) as e:
        samples[:1].plot_distance(metric='jaccard', xlabel='my xlabel', ylabel='my ylabel', title='my title')
    assert 'requires 2 or more' in str(e.value)

    # tooltip with invalid metadata fields or taxids
    with pytest.raises(OneCodexException) as e:
        samples.plot_distance(tooltip='does_not_exist')
    assert 'not found' in str(e.value)


@pytest.mark.parametrize('metric', ('braycurtis', 'jaccard', 'unifrac', 'manhattan'))
def test_plot_mds(ocx, api_data, metric):
    samples = ocx.Samples.where(project='4b53797444f846c4')

    samples.plot_mds(method='pcoa', metric=metric, xlabel='my xlabel', ylabel='my ylabel', title='my title')
    samples.plot_mds(method='smacof', metric=metric, xlabel='my xlabel', ylabel='my ylabel', title='my title')


def test_plot_mds_exceptions(ocx, api_data):
    samples = ocx.Samples.where(project='4b53797444f846c4')

    # expect error if rank is None, since that could lead to weird results
    with pytest.raises(OneCodexException) as e:
        samples.plot_mds(rank=None)
    assert 'specify a rank' in str(e.value)

    # only some metrics allowed
    with pytest.raises(OneCodexException) as e:
        samples.plot_mds(metric='simpson')
    assert 'must be one of' in str(e.value)

    # need more than one analysis
    with pytest.raises(OneCodexException) as e:
        samples[:1].plot_mds(metric='jaccard', xlabel='my xlabel', ylabel='my ylabel', title='my title')
    assert 'requires 2 or more' in str(e.value)

    # tooltip with invalid metadata fields or taxids
    with pytest.raises(OneCodexException) as e:
        samples.plot_mds(tooltip='does_not_exist')
    assert 'not found' in str(e.value)
