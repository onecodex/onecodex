from collections import defaultdict, OrderedDict
from datetime import datetime
import json
import six
import warnings

from onecodex.exceptions import OneCodexException
try:
    from onecodex.analyses import AnalysisMixin
except ImportError:
    class AnalysisMixin(object):
        pass
from onecodex.models import ResourceList


class SampleCollection(ResourceList, AnalysisMixin):
    """A collection of `Samples` or `Classifications` objects with many methods for analysis of
    classifications results.

    Notes
    -----
    Inherits from `ResourceList` to provide a list-like API, and `AnalysisMixin` to provide relevant
    analysis methods.
    """

    def __init__(self, _resource, oc_model, skip_missing=True, label=None, field='auto'):
        self._kwargs = {'skip_missing': skip_missing,
                        'label': label,
                        'field': field}

        super(SampleCollection, self).__init__(_resource, oc_model)

    def _update(self):
        self._cached = {}
        super(SampleCollection, self)._update()

    def _classification_fetch(self, skip_missing=None):
        """Turns a list of objects associated with a classification result into a list of
        Classifications objects.

        Parameters
        ----------
        skip_missing : `bool`
            If an analysis was not successful, exclude it, warn, and keep going

        Returns
        -------
        None, but stores a result in self._cached.
        """
        skip_missing = skip_missing if skip_missing else self._kwargs['skip_missing']

        new_classifications = []

        for a in self._res_list:
            if a.__class__.__name__ == 'Samples':
                c = a.primary_classification
            elif a.__class__.__name__ == 'Classifications':
                c = a
            else:
                raise OneCodexException('{} must be one of: Classifications, Samples')

            if skip_missing and not c.success:
                warnings.warn('Classification {} not successful. Skipping.'.format(c.id))
                continue

            new_classifications.append(c)

        self._cached['classifications'] = new_classifications

    @property
    def primary_classifications(self):
        if 'classifications' not in self._cached:
            self._classification_fetch()

        return self._cached['classifications']

    def _collate_metadata(self, label=None):
        """Turns a list of objects associated with a classification result into a DataFrame of
        metadata.

        Parameters
        ----------
        label : `string` or `callable`
            A metadata field (or function) used to label each analysis. If passing a function, a
            dict containing the metadata for each analysis is passed as the first and only
            positional argument.

        Returns
        -------
        None, but stores a result in self._cached.
        """
        import pandas as pd

        DEFAULT_FIELDS = None
        label = label if label else self._kwargs['label']
        metadata = []

        for c in self.primary_classifications:
            m = c.sample.metadata

            if DEFAULT_FIELDS is None:
                DEFAULT_FIELDS = list(m._resource._schema['properties'].keys())
                DEFAULT_FIELDS.remove('$uri')
                DEFAULT_FIELDS.remove('sample')

            metadatum = {f: getattr(m, f) for f in DEFAULT_FIELDS}
            metadatum['classification_id'] = c.id
            metadatum['sample_id'] = m.sample.id
            metadatum['metadata_id'] = m.id
            metadatum['created_at'] = m.sample.created_at

            if label is None:
                metadatum['_display_name'] = (
                    metadatum['name'] if metadatum['name'] is not None else c.sample.filename
                )
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
                raise NotImplementedError('Must pass a string or function to `label`.')

            metadatum.update(m.custom)
            metadata.append(metadatum)

        if metadata:
            metadata = pd.DataFrame(metadata).set_index('classification_id')

            if all(pd.isnull(metadata['_display_name'])):
                raise OneCodexException('Could not find any labels for `{}`'.format(label))
        else:
            metadata = pd.DataFrame(columns=['classification_id', 'sample_id', 'metadata_id', 'created_at'])

        self._cached['metadata'] = metadata

    @property
    def metadata(self):
        if 'metadata' not in self._cached:
            self._collate_metadata()

        return self._cached['metadata']

    def _collate_results(self, field=None):
        """For a list of objects associated with a classification result, return the results as a
        DataFrame and dict of taxa info.

        Parameters
        ----------
        field : {'readcount_w_children', 'readcount', 'abundance'}
            Which field to use for the abundance/count of a particular taxon in a sample.

            - 'readcount_w_children': total reads of this taxon and all its descendants
            - 'readcount': total reads of this taxon
            - 'abundance': genome size-normalized relative abundances, from shotgun sequencing

        Returns
        -------
        None, but stores a result in self._cached.
        """
        import pandas as pd

        field = field if field else self._kwargs['field']

        if field not in ('auto', 'abundance', 'readcount', 'readcount_w_children'):
            raise OneCodexException('Specified field ({}) not valid.'.format(field))

        # we'll fill these dicts that eventually turn into DataFrames
        df = {
            'classification_id': [c.id for c in self.primary_classifications]
        }

        tax_info = {
            'tax_id': [],
            'name': [],
            'rank': [],
            'parent_tax_id': []
        }

        if field == 'auto':
            field = 'readcount_w_children'

        self._cached['field'] = field

        for c_idx, c in enumerate(self.primary_classifications):
            # pulling results from mainline is the slowest part of the function
            result = c.results()['table']

            # d contains info about a taxon in result, including name, id, counts, rank, etc.
            for d in result:
                d_tax_id = d['tax_id']

                if d_tax_id not in tax_info['tax_id']:
                    for k in ('tax_id', 'name', 'rank', 'parent_tax_id'):
                        tax_info[k].append(d[k])

                    # first time we've seen this taxon, so make a vector for it
                    df[d_tax_id] = [0] * len(self.primary_classifications)

                df[d_tax_id][c_idx] = d[field]

        # format as a Pandas DataFrame
        df = pd.DataFrame(df) \
               .set_index('classification_id') \
               .fillna(0)

        df.columns.name = 'tax_id'

        tax_info = pd.DataFrame(tax_info) \
                     .set_index('tax_id')

        self._cached['results'] = df
        self._cached['taxonomy'] = tax_info

    @property
    def _field(self):
        if 'field' not in self._cached:
            self._collate_results()

        return self._cached['field']

    @property
    def _results(self):
        if 'results' not in self._cached:
            self._collate_results()

        return self._cached['results']

    @property
    def taxonomy(self):
        if 'taxonomy' not in self._cached:
            self._collate_results()

        return self._cached['taxonomy']

    def to_otu(self, biom_id=None):
        """Converts a list of objects associated with a classification result into a `dict` resembling
        an OTU table.

        Parameters
        ----------
        biom_id : `string`, optional
            Optionally specify an `id` field for the generated v1 BIOM file.

        Returns
        -------
        otu_table : `OrderedDict`
            A BIOM OTU table, returned as a Python OrderedDict (can be dumped to JSON)
        """
        otu_format = 'Biological Observation Matrix 1.0.0'

        # Note: This is exact format URL is required by https://github.com/biocore/biom-format
        otu_url = 'http://biom-format.org'

        otu = OrderedDict({'id': biom_id,
                           'format': otu_format,
                           'format_url': otu_url,
                           'type': 'OTU table',
                           'generated_by': 'One Codex API V1',
                           'date': datetime.now().isoformat(),
                           'rows': [],
                           'columns': [],
                           'matrix_type': 'sparse',
                           'matrix_element_type': 'int'})

        rows = defaultdict(dict)

        tax_ids_to_names = {}
        for classification in self.primary_classifications:
            col_id = len(otu['columns'])  # 0 index

            # Re-encoding the JSON is a bit of a hack, but
            # we need a ._to_dict() method that properly
            # resolves references and don't have one at the moment
            columns_entry = {
                'id': str(classification.id),
                'sample_id': str(classification.sample.id),
                'sample_filename': classification.sample.filename,
                'metadata': json.loads(classification.sample.metadata._to_json(include_references=False)),
            }

            otu['columns'].append(columns_entry)
            sample_df = classification.table()

            for row in sample_df.iterrows():
                tax_id = row[1]['tax_id']
                tax_ids_to_names[tax_id] = row[1]['name']
                rows[tax_id][col_id] = int(row[1]['readcount'])

        num_rows = len(rows)
        num_cols = len(otu['columns'])

        otu['shape'] = [num_rows, num_cols]
        otu['data'] = []

        for present_taxa in sorted(rows):
            # add the row entry
            row_id = len(otu['rows'])
            otu['rows'].append({
                'id': present_taxa,
                'metadata': {
                    'taxonomy': tax_ids_to_names[present_taxa],
                }
            })

            for sample_with_hit in rows[present_taxa]:
                counts = rows[present_taxa][sample_with_hit]
                otu['data'].append([row_id, sample_with_hit, counts])

        return otu
