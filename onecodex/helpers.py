import pandas as pd
import six
import warnings

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


def collate_classification_results(classifications, field='readcount_w_children', rank=None,
                                   remove_zeros=True, table_format='wide', multi_index=False,
                                   normalize=False):
    """For a set of classifications, return the results as a Pandas DataFrame and a dict of taxa info.

    Note: The output format is not guaranteed to be stable at this time (i.e., column orderings,
    types, etc. may change).

    classifications (list) -- list of Classifications to collate

    Filtering options:
        remove_zeros (bool) -- remove taxa that are zero in every classification

    Output options:
        table_format ('wide' | 'long')
            - 'wide': rows are classifications, cols are taxa, elements are counts
            - 'long': rows are observations with three cols: classification_id, tax_id, and count

    Options for tabulation of classification results:
        field ('readcount_w_children' | 'readcount' | 'abundance')
            - 'readcount_w_children': total reads of this taxon and all its descendants
            - 'readcount': total reads of this taxon
            - 'abundance': genome size-normalized relative abundances, from shotgun sequencing
        rank (None | 'kingdom' | 'phylum' | 'class' | 'order' | 'family' | 'genus' | 'species')
            - None: include all ranks
            - 'kingdom' or others: restrict analysis to taxa at this rank
        normalize (bool): convert from read counts to relative abundances (each sample sums to 1.0)
    """

    if multi_index:
        warnings.warn('multi_index has been removed--do not use!')

    if field not in ('abundance', 'readcount', 'readcount_w_children'):
        raise OneCodexException('Specified field ({}) not valid.'.format(field))

    classifications = [c for c in classifications if c.success is True]

    # we'll fill a dict that eventually turn into a pandas df
    df = {
        'classification_id': [c.id for c in classifications]
    }

    tax_info = {}

    for c_idx, c in enumerate(classifications):
        # pulling results from mainline is the slowest part of the function
        result = c.results()['table']

        # d contains info about a taxon in result, including name, id, counts, rank, etc.
        for d in result:
            d_tax_id = d['tax_id']

            if d_tax_id not in tax_info:
                # only process this taxon if it's the correct rank
                if rank is not None and d['rank'] != rank:
                    continue

                tax_info[d_tax_id] = {k: d[k] for k in ('name', 'rank', 'parent_tax_id')}

                # first time we've seen this taxon, so make a vector for it
                df[d_tax_id] = [0] * len(classifications)

            df[d_tax_id][c_idx] = d[field]

    # no results collated--was an incorrect rank specified?
    if len(tax_info) == 0:
        raise OneCodexException('No results collated. Is taxonomic rank ({}) valid?'.format(rank))

    # format as a Pandas DataFrame
    df = pd.DataFrame(df) \
           .set_index('classification_id') \
           .fillna(0)

    # normalize
    if normalize:
        df = df.div(df.sum(axis=1), axis=0)

    # remove columns (tax_ids) with no values that are > 0
    if remove_zeros:
        df = df.loc[:, (df != 0).any(axis=0)]

    if table_format == 'long':
        long_df = {
            'classification_id': [],
            'tax_id': [],
            field: []
        }

        for t_id in df:
            for c_id, count in df[t_id].iteritems():
                long_df['classification_id'].append(c_id)
                long_df['tax_id'].append(t_id)
                long_df[field].append(count)

        df = pd.DataFrame(long_df)

    return df, tax_info
