from functools import partial
import pandas as pd

import six

from onecodex.exceptions import OneCodexException
from onecodex.models import Metadata
from onecodex.viz import VizPCAMixin, VizHeatmapMixin, VizMetadataMixin, VizDistanceMixin


# force persistence of our additional taxonomy and metadata dataframe properties
class ResultsDataFrame(pd.DataFrame):
    _metadata = ['ocx_rank', 'ocx_field', 'ocx_taxonomy', 'ocx_metadata']

    def __init__(self, data=None, index=None, columns=None, dtype=None, copy=False, ocx_data={}):
        self.ocx_rank = ocx_data.get('ocx_rank', None)
        self.ocx_field = ocx_data.get('ocx_field', None)
        self.ocx_taxonomy = ocx_data.get('ocx_taxonomy', None)
        self.ocx_metadata = ocx_data.get('ocx_metadata', None)

        pd.DataFrame.__init__(self, data=data, index=index, columns=columns, dtype=dtype, copy=copy)

    @property
    def _constructor(self):
        # we explicitly do *not* pass rank on to manipulated ResultsDataFrames. we don't know
        # how the data has been manipulated, and it may no longer be accurate
        ocx_data = {
            'ocx_rank': None,
            'ocx_field': self.ocx_field,
            'ocx_taxonomy': self.ocx_taxonomy,
            'ocx_metadata': self.ocx_metadata
        }

        return partial(ResultsDataFrame, ocx_data=ocx_data)

    @property
    def _constructor_sliced(self):
        # we explicitly do *not* pass rank on to manipulated ResultsDataFrames. we don't know
        # how the data has been manipulated, and it may no longer be accurate
        ocx_data = {
            'ocx_rank': None,
            'ocx_field': self.ocx_field,
            'ocx_taxonomy': self.ocx_taxonomy,
            'ocx_metadata': self.ocx_metadata
        }

        return partial(ResultsSeries, ocx_data=ocx_data)


class ResultsSeries(pd.Series):
    _metadata = ['name', 'ocx_rank', 'ocx_field', 'ocx_taxonomy', 'ocx_metadata']

    def __init__(self, data=None, index=None, dtype=None, name=None, copy=False, fastpath=False, ocx_data={}):
        self.ocx_rank = ocx_data.get('ocx_rank', None)
        self.ocx_field = ocx_data.get('ocx_field', None)
        self.ocx_taxonomy = ocx_data.get('ocx_taxonomy', None)
        self.ocx_metadata = ocx_data.get('ocx_metadata', None)

        pd.Series.__init__(self, data=data, index=index, dtype=dtype, name=name, copy=copy, fastpath=fastpath)

    @property
    def _constructor(self):
        # we explicitly do *not* pass rank on to manipulated ResultsDataFrames. we don't know
        # how the data has been manipulated, and it may no longer be accurate
        ocx_data = {
            'ocx_rank': None,
            'ocx_field': self.ocx_field,
            'ocx_taxonomy': self.ocx_taxonomy,
            'ocx_metadata': self.ocx_metadata
        }

        return partial(ResultsSeries, ocx_data=ocx_data)

    @property
    def _constructor_expanddim(self):
        # we explicitly do *not* pass rank on to manipulated ResultsDataFrames. we don't know
        # how the data has been manipulated, and it may no longer be accurate
        ocx_data = {
            'ocx_rank': None,
            'ocx_field': self.ocx_field,
            'ocx_taxonomy': self.ocx_taxonomy,
            'ocx_metadata': self.ocx_metadata
        }

        return partial(ResultsDataFrame, ocx_data=ocx_data)


class AnalysisMethods(VizPCAMixin, VizHeatmapMixin, VizMetadataMixin, VizDistanceMixin):
    _immutable = ['field', '_results', '_metadata', '_taxonomy']

    def __init__(self, results, metadata, taxonomy, field):
        self.__setattr__('_results', results, force=True)
        self.__setattr__('_metadata', metadata, force=True)
        self.__setattr__('_taxonomy', taxonomy, force=True)
        self.__setattr__('field', field, force=True)

    def __setattr__(self, name, value, force=False):
        if name in self._immutable and not force:
            raise OneCodexException('{} is an immutable property and can not be changed'.format(name))

        super(AnalysisMethods, self).__setattr__(name, value)

    def _get_auto_rank(self, rank):
        """Tries to figure out what rank we should use for analyses, mainly called by results()"""

        if rank == 'auto':
            # if we're an accessor for a ResultsDataFrame, use its ocx_rank property
            if isinstance(self, OneCodexAccessor):
                return self._results.ocx_rank

            if self.field == 'abundance':
                return 'species'
            else:
                return 'genus'
        else:
            return rank

    def _guess_normalized(self):
        # it's possible that the _results df has already been normalized, which can cause some
        # methods to fail. we must guess whether this is the case

        return bool((self._results.sum(axis=1).round(4) == 1.0).all())

    def magic_metadata_fetch(self, metadata_fields):
        """Takes a list of metadata fields, some of which can contain taxon names or taxon IDs, and
        returns a DataFrame with magically transformed data that can be used for plotting.
        """

        help_metadata = ', '.join(self._metadata.keys())
        magic_metadata = pd.DataFrame({'classification_id': self._results.index}) \
                           .set_index('classification_id')

        # if we magically rename fields, keep track
        magic_fields = {}

        for f in set([f for f in metadata_fields if f]):
            if isinstance(f, tuple):
                # joined categorical metadata
                for field in f:
                    if field not in self._metadata:
                        raise OneCodexException(
                            'Field {} not found. Choose from: {}'.format(field, help_metadata)
                        )

                    if not (pd.api.types.is_bool_dtype(self._metadata[field]) or  # noqa
                            pd.api.types.is_categorical_dtype(self._metadata[field]) or  # noqa
                            pd.api.types.is_object_dtype(self._metadata[field])):
                        raise OneCodexException(
                            'When specifying multiple metadata fields, all must be categorical'
                        )

                # concatenate the columns together with underscores
                composite_field = '_'.join(f)
                magic_metadata[composite_field] = ''
                magic_metadata[composite_field] = magic_metadata[composite_field].str.cat(
                    [self._metadata[field].astype(str) for field in f], sep='_'
                ).str.lstrip('_')
                magic_fields[f] = composite_field
            else:
                str_f = str(f)

                if str_f == 'Label':
                    # it's our magic keyword that means _display_name
                    magic_metadata[f] = self._metadata['_display_name']
                    magic_fields[f] = str_f
                elif str_f in self._metadata:
                    # exactly matches existing metadata field
                    magic_metadata[f] = self._metadata[str_f]
                    magic_fields[f] = str_f
                elif str_f in self._results.keys():
                    # is a tax_id
                    tax_name = self._taxonomy['name'][str_f]

                    # report within-rank abundance
                    df = self.results(rank=self._taxonomy['rank'][str_f])

                    renamed_field = '{} ({})'.format(tax_name, str_f)
                    magic_metadata[renamed_field] = df[str_f]
                    magic_fields[f] = renamed_field
                else:
                    # try to match it up with a taxon name
                    hits = []

                    # don't both searching if the query is really short
                    if len(str_f) > 4:
                        for tax_id, tax_name in zip(self._taxonomy.index, self._taxonomy['name']):
                            # if it's an exact match, use that and skip the rest
                            if str_f.lower() == tax_name.lower():
                                hits = [(tax_id, tax_name)]
                                break
                            # otherwise, keep trying to match
                            elif str_f.lower() in tax_name.lower():
                                hits.append((tax_id, tax_name))

                        # take the hit with the lowest tax_id
                        hits = sorted(hits, key=lambda x: int(x[0]))

                    if hits:
                        # report within-rank abundance
                        df = self.results(rank=self._taxonomy['rank'][hits[0][0]])

                        renamed_field = '{} ({})'.format(hits[0][1], hits[0][0])
                        magic_metadata[renamed_field] = df[hits[0][0]]
                        magic_fields[f] = renamed_field
                    else:
                        # matched nothing
                        raise OneCodexException(
                            'Field or taxon {} not found. Choose from: {}'.format(str_f, help_metadata)
                        )

        return magic_metadata, magic_fields

    def results(self, rank='auto', top_n=None, threshold=None,
                remove_zeros=True, normalize='auto',
                table_format='wide', **kwargs):
        """
        Filtering options (performed consecutively in the order given below):
            rank (None | 'auto' | 'kingdom' | 'phylum' | 'class' | 'order' | 'family' | 'genus' | 'species')
                - None: include all ranks
                - 'auto': choose automatically based on fields
                - 'kingdom' or others: only return taxa at this rank
            remove_zeros (bool) -- remove taxa that are zero in every classification
            threshold (float) -- only return taxa more abundant than this threshold
            top_n (int) -- return the top N most abundant taxa

        Output options:
            table_format ('wide' | 'long')
                - 'wide': rows are classifications, cols are taxa, elements are counts
                - 'long': rows are observations with 3 cols: classification_id, tax_id, and count
            normalize (bool): convert read counts to relative abundances (each sample sums to 1.0)
        """

        rank = self._get_auto_rank(rank)
        df = self._results.copy()

        # subset by taxa
        if rank:
            tax_ids_to_keep = []

            for tax_id in df.keys():
                if self._taxonomy['rank'][tax_id] == rank:
                    tax_ids_to_keep.append(tax_id)

            if len(tax_ids_to_keep) == 0:
                raise OneCodexException('No taxa kept--is rank ({}) correct?'.format(rank))

            df = df.loc[:, tax_ids_to_keep]

        # normalize
        if normalize is False and self._guess_normalized():
            raise OneCodexException('Data has already been normalized and this can not be undone.')

        if normalize is True or (normalize == 'auto' and rank is not None and self.field != 'abundance'):
            df = df.div(df.sum(axis=1), axis=0)

        # remove columns (tax_ids) with no values that are > 0
        if remove_zeros:
            df = df.loc[:, (df != 0).any(axis=0)]

        # restrict to taxa appearing in one or more samples at the given threshold
        if threshold:
            df = df.loc[:, df.max() >= threshold]

        # restrict to N most abundant taxa
        if top_n:
            idx = df.sum(axis=0).sort_values(ascending=False).head(top_n).index
            df = df.loc[:, idx]

        # generate long-format table
        if table_format == 'long':
            long_df = {
                'classification_id': [],
                'tax_id': [],
                self.field: []
            }

            for t_id in df:
                for c_id, count in df[t_id].iteritems():
                    long_df['classification_id'].append(c_id)
                    long_df['tax_id'].append(t_id)
                    long_df[self.field].append(count)

            results_df = ResultsDataFrame(long_df)
        elif table_format == 'wide':
            results_df = ResultsDataFrame(df)
        else:
            raise OneCodexException('table_format must be one of: long, wide')

        results_df.ocx_rank = rank
        results_df.ocx_field = self.field
        results_df.ocx_taxonomy = self.taxonomy()
        results_df.ocx_metadata = self.metadata()

        return results_df

    def metadata(self):
        return self._metadata.copy()

    def taxonomy(self):
        return self._taxonomy.copy()


class EnhancedSampleCollection(AnalysisMethods):
    def _update(self):
        # whenever a change is made to the samples in this list (i.e., removing some of them), we must update
        # the metadata and results dataframes we're saving in this class. use the kwargs that were originally
        # passed when this object was instantiated.
        super(EnhancedSampleCollection, self)._update()

        if self._resource:
            self.collate_metadata(**{'label': self._kwargs['label']} if 'label' in self._kwargs else {})
            self.collate_results(**{'field': self._kwargs['field']} if 'field' in self._kwargs else {})

    def collate_metadata(self, label=None):
        """Turns a list of objects associated with a classification result into a DataFrame of
        metadata.

        analyses (list) -- list of Samples, Classifications, or Analyses objects
        label (string | function) -- metadata field (or function) used to label each analysis. if
            passing a function, a dict containing the metadata for each analysis is passed as the
            first and only positional argument.
        """

        metadata = []

        DEFAULT_FIELDS = list(Metadata._resource._schema['properties'].keys())
        DEFAULT_FIELDS.remove('$uri')
        DEFAULT_FIELDS.remove('sample')

        for c in self._classifications:
            m = c.sample.metadata

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

        metadata = pd.DataFrame(metadata).set_index('classification_id')

        if all(pd.isnull(metadata['_display_name'])):
            raise OneCodexException('Could not find any labels for `{}`'.format(label))

        self.__setattr__('_metadata', metadata, force=True)

    def collate_results(self, field='auto'):
        """For a list of objects associated with a classification result, return the results as a
        DataFrame and dict of taxa info.

        field ('readcount_w_children' | 'readcount' | 'abundance')
            - 'readcount_w_children': total reads of this taxon and all its descendants
            - 'readcount': total reads of this taxon
            - 'abundance': genome size-normalized relative abundances, from shotgun sequencing
        """

        if field not in ('auto', 'abundance', 'readcount', 'readcount_w_children'):
            raise OneCodexException('Specified field ({}) not valid.'.format(field))

        # we'll fill these dicts that eventually turn into DataFrames
        df = {
            'classification_id': [c.id for c in self._classifications]
        }

        tax_info = {
            'tax_id': [],
            'name': [],
            'rank': [],
            'parent_tax_id': []
        }

        if field == 'auto':
            field = 'readcount_w_children'

        self.__setattr__('field', field, force=True)

        for c_idx, c in enumerate(self._classifications):
            # pulling results from mainline is the slowest part of the function
            result = c.results()['table']

            # d contains info about a taxon in result, including name, id, counts, rank, etc.
            for d in result:
                d_tax_id = d['tax_id']

                if d_tax_id not in tax_info['tax_id']:
                    for k in ('tax_id', 'name', 'rank', 'parent_tax_id'):
                        tax_info[k].append(d[k])

                    # first time we've seen this taxon, so make a vector for it
                    df[d_tax_id] = [0] * len(self._classifications)

                df[d_tax_id][c_idx] = d[field]

        # format as a Pandas DataFrame
        df = pd.DataFrame(df) \
               .set_index('classification_id') \
               .fillna(0)

        tax_info = pd.DataFrame(tax_info) \
                     .set_index('tax_id')

        self.__setattr__('_results', df, force=True)
        self.__setattr__('_taxonomy', tax_info, force=True)


@pd.api.extensions.register_dataframe_accessor('ocx')
class OneCodexAccessor(AnalysisMethods):
    def __init__(self, pandas_obj):
        # copy data from the ResultsDataFrame to a new instance of AnalysisMethods
        super(OneCodexAccessor, self).__init__(
            pandas_obj,
            pandas_obj.ocx_metadata,
            pandas_obj.ocx_taxonomy,
            pandas_obj.ocx_field
        )

        # prune back _taxonomy df to contain only taxa present in the ResultsDataFrame (and parents)
        tree = self.tree_build()
        tree = self.tree_prune_tax_ids(tree, self._results.keys())

        tax_ids_to_keep = [x.name for x in tree.traverse()]

        super(OneCodexAccessor, self).__setattr__('_taxonomy', self._taxonomy.loc[tax_ids_to_keep], force=True)

        # similarly restrict _metadata df to contain only data relevant to samples currently in ResultsDataFrame
        super(OneCodexAccessor, self).__setattr__('_metadata', self._metadata.loc[self._results.index], force=True)
