import numpy as np
import pandas as pd

from onecodex.exceptions import OneCodexException
from onecodex.models import Analyses, Classifications, Samples


def normalize_analyses(analyses, skip_missing=True):
    normed_analyses = []
    metadata = []

    for a in analyses:
        if isinstance(a, Samples):
            c = a.primary_classification
            m = a.metadata
        elif isinstance(a, Classifications):
            c = a
            m = a.sample.metadata
        elif isinstance(a, Analyses):
            if a.analysis_type != 'classification':
                raise OneCodexException('{} is not a classification'.format(a.id))
            c = Classifications(a._resource._client.Classifications.fetch(a.id))
            m = a.sample.metadata

        if skip_missing and not c.success:
            continue

        normed_analyses.append(c)
        metadata.append(m)

    return normed_analyses, metadata


def collate_analysis_results(analyses, field='abundance'):
    """For a set of analyses, return the results as a Pandas DataFrame."""
    assert field in ['abundance', 'readcount', 'readcount_w_children']

    # Keep track of all of the microbial abundances
    dat = []
    titles = []

    # Keep track of information for each tax_id
    tax_id_info = {}

    # Get results for each of the Sample objects that are passed in
    for a in analyses:
        if a.success is False:
            dat.append({})
            titles.append(a.id)
            continue

        # Get the results in table format
        result = a.results()['table']

        # Record the information (name, parent) for each organism by  its tax ID
        for d in result:
            if d['tax_id'] not in tax_id_info:
                tax_id_info[d['tax_id']] = {k: d[k] for k in ['name', 'rank', 'parent_tax_id']}

        # Reformat detection infromation as dict of {taxid: value}
        result = {d['tax_id']: d[field] for d in result if field in d}

        # Remove entries without the specified field
        result = {taxid: value for taxid, value in result.items() if value is not None}

        # Save the set of microbial abundances
        dat.append(result)
        titles.append(a.id)

    if len(dat) == 0:
        return None

    # Format as a Pandas DataFrame
    df = pd.DataFrame(dat)
    df.index = titles

    # fill in missing values
    df = df.T.fillna(0)

    # add an index with the tax ids name
    df.index.name = 'tax_id'
    df['tax_name'] = df.index.map(lambda tid: tax_id_info[tid]['name'])
    df['tax_rank'] = df.index.map(lambda tid: tax_id_info[tid]['rank'])
    df.set_index(['tax_name', 'tax_rank'], inplace=True, append=True)

    # Remove columns (tax_ids) with no values that are > 0
    df = df.T.loc[:, df.T.sum() > 0]

    return df
