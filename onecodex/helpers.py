from functools import partial
import pandas as pd

from onecodex.exceptions import OneCodexException
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
    def __init__(self, results, metadata, taxonomy, field):
        self._cached = {'results': results,
                        'metadata': metadata,
                        'taxonomy': taxonomy,
                        'field': field}

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

    def _metadata_fetch(self, metadata_fields):
        """Takes a list of metadata fields, some of which can contain taxon names or taxon IDs, and
        returns a DataFrame with magically transformed data that can be used for plotting.
        """
        help_metadata = ', '.join(self.metadata.keys())
        magic_metadata = pd.DataFrame({'classification_id': self._results.index}) \
                           .set_index('classification_id')

        # if we magically rename fields, keep track
        magic_fields = {}

        for f in set([f for f in metadata_fields if f]):
            if isinstance(f, tuple):
                # joined categorical metadata
                for field in f:
                    if field not in self.metadata:
                        raise OneCodexException(
                            'Field {} not found. Choose from: {}'.format(field, help_metadata)
                        )

                    if not (pd.api.types.is_bool_dtype(self.metadata[field]) or  # noqa
                            pd.api.types.is_categorical_dtype(self.metadata[field]) or  # noqa
                            pd.api.types.is_object_dtype(self.metadata[field])):
                        raise OneCodexException(
                            'When specifying multiple metadata fields, all must be categorical'
                        )

                # concatenate the columns together with underscores
                composite_field = '_'.join(f)
                magic_metadata[composite_field] = ''
                magic_metadata[composite_field] = magic_metadata[composite_field].str.cat(
                    [self.metadata[field].astype(str) for field in f], sep='_'
                ).str.lstrip('_')
                magic_fields[f] = composite_field
            else:
                str_f = str(f)

                if str_f == 'Label':
                    # it's our magic keyword that means _display_name
                    magic_metadata[f] = self.metadata['_display_name']
                    magic_fields[f] = str_f
                elif str_f in self.metadata:
                    # exactly matches existing metadata field
                    magic_metadata[f] = self.metadata[str_f]
                    magic_fields[f] = str_f
                elif str_f in self._results.keys():
                    # is a tax_id
                    tax_name = self.taxonomy['name'][str_f]

                    # report within-rank abundance
                    df = self.results(rank=self.taxonomy['rank'][str_f])

                    renamed_field = '{} ({})'.format(tax_name, str_f)
                    magic_metadata[renamed_field] = df[str_f]
                    magic_fields[f] = renamed_field
                else:
                    # try to match it up with a taxon name
                    hits = []

                    # don't both searching if the query is really short
                    if len(str_f) > 4:
                        for tax_id, tax_name in zip(self.taxonomy.index, self.taxonomy['name']):
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
                        df = self.results(rank=self.taxonomy['rank'][hits[0][0]])

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
                table_format='wide'):
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
                if self.taxonomy['rank'][tax_id] == rank:
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
        results_df.ocx_taxonomy = self.taxonomy.copy()
        results_df.ocx_metadata = self.metadata.copy()

        return results_df

    @property
    def field(self):
        return self._cached['field']

    @property
    def metadata(self):
        return self._cached['metadata']

    @property
    def _results(self):
        return self._cached['results']

    @property
    def taxonomy(self):
        return self._cached['taxonomy']


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

        self._cached['taxonomy'] = self.taxonomy.loc[tax_ids_to_keep]

        # similarly restrict _metadata df to contain only data relevant to samples currently in ResultsDataFrame
        self._cached['metadata'] = self.metadata.loc[self._results.index]
