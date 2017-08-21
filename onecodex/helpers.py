import pandas as pd
import warnings

import six

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
                                   rank=None, remove_zeros=True):
    """For a set of classifications, return the results as a Pandas DataFrame."""
    assert field in ['abundance', 'readcount', 'readcount_w_children']

    # Keep track of all of the microbial abundances
    dat = []
    titles = []

    # Keep track of information for each tax_id
    tax_id_info = {}

    # Get results for each of the Sample objects that are passed in
    # Note: Warnings are handled in the `normalize_classifications` function
    for c in classifications:
        if c.success is False:
            dat.append({})
            titles.append(c.id)
            continue

        # Get the results in table format
        result = c.results()['table']

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
        titles.append(c.id)

    if len(dat) == 0:
        return None

    # Format as a Pandas DataFrame
    # TODO: Optimize this; unexpectedly slow
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
    if remove_zeros:
        df = df.T.loc[:, df.T.sum() > 0]
    else:
        df = df.T

    # TODO: Move this up further to optimize for smaller dataframe above
    if rank is not None:
        df = df.loc[:, [i[2] == rank for i in df.columns]]

    return df
