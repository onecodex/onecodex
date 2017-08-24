import numpy as np
import pandas as pd
import warnings

import six

from collections import OrderedDict

from onecodex.exceptions import OneCodexException
from onecodex.models import Analyses, Classifications, Samples, Metadata


def normalize_classifications(analyses, label=None, skip_missing=True, warn=True):
    """Get a list of classifications and metadata from input `analyses`.

    These "analyses" may be Samples, Analyses, Classifications, or a mix
    thereof. Metadata is taken off of the Classification.sample.metadata object.
    """
    normed_classifications = []
    metadata = []

    DEFAULT_FIELDS = list(Metadata._resource._schema['properties'].keys())
    DEFAULT_FIELDS.remove('$uri')
    DEFAULT_FIELDS.remove('sample')

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
            if warn:
                warnings.warn('Classification {} not successful. Skipping.'.format(c.id))
            continue

        normed_classifications.append(c)
        metadatum = {f: getattr(m, f) for f in DEFAULT_FIELDS}

        # Add sample_id, metadata_id, and sample created_at
        metadatum['sample_id'] = m.sample.id
        metadatum['metadata_id'] = m.id
        metadatum['created_at'] = m.sample.created_at

        if label is None:
            metadatum['_display_name'] = (metadatum['name'] if metadatum['name'] is not None
                                          else c.sample.filename)
        elif isinstance(label, six.string_types):
            if label in metadatum:
                metadatum['_display_name'] = metadatum[label]
            elif label in m.custom:
                metadatum['_display_name'] = m.custom[label]
            else:
                metadatum['_display_name'] = None
        elif callable(label):
            metadatum['_display_name'] = label(m)
        else:
            raise NotImplementedError('Must pass a string or lambda function to `label`.')

        metadatum.update(m.custom)
        metadata.append(metadatum)

    metadata = pd.DataFrame(metadata)
    if all(pd.isnull(metadata['_display_name'])):
        raise OneCodexException('Could not find any labels for `{}`'.format(label))

    return normed_classifications, metadata


def collate_classification_results(classifications, field='readcount_w_children',
                                   rank=None, remove_zeros=True, multi_index=False,
                                   normalize=False):
    """For a set of classifications, return the results as a Pandas DataFrame and a dict of taxa info.

    Note: The output format is not guaranteed to be stable at this time (i.e.,
    column orderings, types, etc. may change).
    """
    assert field in ['abundance', 'readcount', 'readcount_w_children']
    TYPE_MAPPING = {
        'readcount': np.int64,
        'readcount_w_children': np.int64,
        'abundance': np.float64,
    }

    # Keep track of all of the microbial abundances
    titles = []
    tax_id_info = OrderedDict()
    column_ix = 0
    for c in classifications:
        if c.success is False:
            continue

        # Get the results in table format
        result = c.results()['table']

        # Record the information (name, parent) for each organism by  its tax ID
        for d in result:
            if d['tax_id'] not in tax_id_info:
                tax_id_info[d['tax_id']] = {k: d[k] for k in ['name', 'rank', 'parent_tax_id']}
                tax_id_info[d['tax_id']]['column'] = column_ix
                column_ix += 1

    arr = np.zeros((len(classifications), len(tax_id_info)), dtype=TYPE_MAPPING[field])
    for ix, c in enumerate(classifications):
        titles.append(c.id)
        result = c.results()['table']
        for d in result:
            arr[ix, tax_id_info[d['tax_id']]['column']] = d[field]

    # Format as a Pandas DataFrame
    df = pd.DataFrame(arr, columns=tax_id_info.keys())
    df.index = titles

    # fill in missing values
    df = df.T.fillna(0)

    # add an index with the tax ids name
    df.index.name = 'tax_id'
    if multi_index:
        df['tax_name'] = df.index.map(lambda tid: tax_id_info[tid]['name'])
        df['tax_rank'] = df.index.map(lambda tid: tax_id_info[tid]['rank'])
        df.set_index(['tax_name', 'tax_rank'], inplace=True, append=True)
        df.sort_index(inplace=True)

    # Subset to rank as appropriate
    if rank is not None:
        df = df.loc[[k for k, v in tax_id_info.items() if v['rank'] == rank], :]

    # Normalize
    if normalize:
        df = df.div(df.sum(axis=1), axis=0)

    # Remove columns (tax_ids) with no values that are > 0
    if remove_zeros:
        df = df.T.loc[:, df.T.sum() > 0]
    else:
        df = df.T

    # Set transposed index name
    df.index.name = 'classification_id'

    return df, tax_id_info
