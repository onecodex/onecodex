from onecodex.helpers import collate_analysis_results


def test_result_collation(ocx, api_data):
    analyses = [ocx.Classifications.get('f9e4a5506b154953')]
    results = collate_analysis_results(analyses)

    assert 'f9e4a5506b154953' in results.index
    assert len(results.loc['f9e4a5506b154953']) > 0
    assert 'tax_id' in results.columns.names and 'tax_name' in results.columns.names
