import pytest
import skbio

from onecodex.exceptions import OneCodexException
from onecodex.helpers import collate_classification_results
from onecodex.viz import plot_distance, plot_heatmap, plot_metadata, plot_pca


def test_result_collation(ocx, api_data):
    analyses = [ocx.Classifications.get('45a573fb7833449a')]
    results, tax_info = collate_classification_results(analyses)

    assert '45a573fb7833449a' in results.index
    assert len(results.loc['45a573fb7833449a']) > 0
    assert 'tax_id' in results.columns.names


# Note: Need a better plotting setup, these tests just ensure
# we don't blow up / raise any exceptions in the plotting
def test_plot_metadata(ocx, api_data):
    analyses = [ocx.Classifications.get('45a573fb7833449a')]
    samples = [a.sample for a in analyses]

    plot_metadata(analyses, statistic='simpson')

    # Should resolve via samples OK too; date objects should get coerced OK
    plot_metadata(samples, metadata='date_sequenced', statistic='simpson')

    # And bools
    plot_metadata(samples, metadata='starred', statistic='simpson')

    # Numbers too
    samples[0].metadata.custom['my_int'] = 10
    plot_metadata(samples, metadata='my_int', statistic='simpson',
                  xlabel='X Axis', ylabel='Y Axis', title='Title')


def test_plot_metadata_exceptions(ocx, api_data):
    analyses = [ocx.Classifications.get('45a573fb7833449a')]

    # Must specify a statistic or tax_id
    with pytest.raises(OneCodexException):
        plot_metadata(analyses)

    # Cannot specify both
    with pytest.raises(OneCodexException):
        plot_metadata(analyses, statistic='simpson', tax_id=562)

    # Metadata field must exist
    with pytest.raises(OneCodexException):
        plot_metadata(analyses, metadata='invalid_field')

    # Statistic must be valid
    with pytest.raises(NotImplementedError):
        plot_metadata(analyses, statistic='simpson2')


def test_plot_pca(ocx, api_data):
    # Requires >1 classification or blows up
    analyses = [ocx.Classifications.get('45a573fb7833449a'),
                ocx.Classifications.get('593601a797914cbf')]

    analyses = analyses + analyses  # Note a uniqueness check may eventually render this an error
    plot_pca(analyses, hue='platform', org_vectors=1)


def test_plot_pca_exceptions(ocx, api_data):
    # Requires >1 classification or blows up
    analyses = [ocx.Classifications.get('45a573fb7833449a')]
    with pytest.raises(OneCodexException):
        plot_pca(analyses, hue='platform')


def test_plot_heatmap(ocx, api_data):
    analyses = [ocx.Classifications.get('45a573fb7833449a'),
                ocx.Classifications.get('593601a797914cbf')]

    plot_heatmap(analyses, top_n=10)
    plot_heatmap(analyses, threshold=0.1)
    plot_heatmap(analyses, top_n=10, threshold=0.1)
    plot_heatmap(analyses, top_n=10, threshold=0.1, xlabel='X Axis', ylabel='Y Axis', title='Title')


def test_plot_heatmap_exceptions(ocx, api_data):
    analyses = [ocx.Classifications.get('45a573fb7833449a')]

    # Requires multiple analyses
    with pytest.raises(OneCodexException):
        plot_heatmap(analyses)


def test_plot_distances(ocx, api_data):
    analyses = [ocx.Classifications.get('45a573fb7833449a'),
                ocx.Classifications.get('593601a797914cbf')]

    # Invalid distance
    with pytest.raises(OneCodexException):
        plot_distance(analyses, metric='simpson')

    plot_distance(analyses, metric='jaccard', xlabel='X Axis', ylabel='Y Axis', title='Title')


def test_plot_distances_exceptions(ocx, api_data):
    analyses = [ocx.Classifications.get('45a573fb7833449a')]

    with pytest.raises(OneCodexException):
        plot_distance(analyses, metric='jaccard', xlabel='X Axis', ylabel='Y Axis', title='Title')

    # Can't have just duplicate data, blows up skbio
    analyses = analyses + analyses
    with pytest.raises(skbio.stats.distance.DissimilarityMatrixError):
        plot_distance(analyses, metric='jaccard', xlabel='X Axis', ylabel='Y Axis', title='Title')
