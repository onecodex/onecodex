from functools import partial
import pandas as pd

from onecodex.exceptions import OneCodexException
from onecodex.viz import VizPCAMixin, VizHeatmapMixin, VizMetadataMixin, VizDistanceMixin


# force persistence of our additional taxonomy and metadata dataframe properties
class ClassificationsDataFrame(pd.DataFrame):
    _metadata = ['ocx_rank', 'ocx_field', 'ocx_taxonomy', 'ocx_metadata']

    def __init__(self, data=None, index=None, columns=None, dtype=None, copy=False, ocx_rank=None,
                 ocx_field=None, ocx_taxonomy=None, ocx_metadata=None):
        self.ocx_rank = ocx_rank
        self.ocx_field = ocx_field
        self.ocx_taxonomy = ocx_taxonomy
        self.ocx_metadata = ocx_metadata

        pd.DataFrame.__init__(self, data=data, index=index, columns=columns, dtype=dtype, copy=copy)

    @property
    def _constructor(self):
        # we explicitly do *not* pass rank on to manipulated ClassificationsDataFrame. we don't know
        # how the data has been manipulated, and it may no longer be accurate
        return partial(ClassificationsDataFrame, ocx_rank=None, ocx_field=self.ocx_field,
                       ocx_taxonomy=self.ocx_taxonomy, ocx_metadata=self.ocx_metadata)

    @property
    def _constructor_sliced(self):
        # we explicitly do *not* pass rank on to manipulated ClassificationsDataFrame. we don't know
        # how the data has been manipulated, and it may no longer be accurate
        return partial(ClassificationsSeries, ocx_rank=None, ocx_field=self.ocx_field,
                       ocx_taxonomy=self.ocx_taxonomy, ocx_metadata=self.ocx_metadata)


class ClassificationsSeries(pd.Series):
    # 'name' is a piece of metadata specified by pd.Series--it's not ours
    _metadata = ['name', 'ocx_rank', 'ocx_field', 'ocx_taxonomy', 'ocx_metadata']

    def __init__(self, data=None, index=None, dtype=None, name=None, copy=False, fastpath=False,
                 ocx_rank=None, ocx_field=None, ocx_taxonomy=None, ocx_metadata=None):
        self.ocx_rank = ocx_rank
        self.ocx_field = ocx_field
        self.ocx_taxonomy = ocx_taxonomy
        self.ocx_metadata = ocx_metadata

        pd.Series.__init__(self, data=data, index=index, dtype=dtype, name=name, copy=copy, fastpath=fastpath)

    @property
    def _constructor(self):
        # we explicitly do *not* pass rank on to manipulated ClassificationsDataFrames. we don't know
        # how the data has been manipulated, and it may no longer be accurate
        return partial(ClassificationsSeries, ocx_rank=None, ocx_field=self.ocx_field,
                       ocx_taxonomy=self.ocx_taxonomy, ocx_metadata=self.ocx_metadata)

    @property
    def _constructor_expanddim(self):
        # we explicitly do *not* pass rank on to manipulated ClassificationsDataFrame. we don't know
        # how the data has been manipulated, and it may no longer be accurate
        return partial(ClassificationsDataFrame, ocx_rank=None, ocx_field=self.ocx_field,
                       ocx_taxonomy=self.ocx_taxonomy, ocx_metadata=self.ocx_metadata)


class AnalysisMixin(VizPCAMixin, VizHeatmapMixin, VizMetadataMixin, VizDistanceMixin):
    def _get_auto_rank(self, rank):
        """Tries to figure out what rank we should use for analyses, mainly called by results()"""

        if rank == 'auto':
            # if we're an accessor for a ClassificationsDataFrame, use its _rank property
            if isinstance(self, OneCodexAccessor):
                return self._rank

            if self._field == 'abundance':
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
                    df = self.to_df(rank=self.taxonomy['rank'][str_f])

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
                        df = self.to_df(rank=self.taxonomy['rank'][hits[0][0]])

                        renamed_field = '{} ({})'.format(hits[0][1], hits[0][0])
                        magic_metadata[renamed_field] = df[hits[0][0]]
                        magic_fields[f] = renamed_field
                    else:
                        # matched nothing
                        raise OneCodexException(
                            'Field or taxon {} not found. Choose from: {}'.format(str_f, help_metadata)
                        )

        return magic_metadata, magic_fields

    def to_df(self, rank='auto', top_n=None, threshold=None, remove_zeros=True, normalize='auto',
              table_format='wide'):
        """Takes the ClassificationsDataFrame associated with these samples, or SampleCollection, does some
        filtering, and returns a ClassificationsDataFrame.

        Parameters
        ----------
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.
        top_n : `integer`, optional
            Return only the top N most abundant taxa.
        threshold : `float`, optional
            Return only taxa more abundant than this threshold in one or more samples.
        remove_zeros : `bool`, optional
            Do not return taxa that have zero abundance in every sample.
        normalize : {'auto', True, False}
            Convert read counts to relative abundances (each sample sums to 1.0).
        table_format : {'long', 'wide'}
            If wide, rows are classifications, cols are taxa, elements are counts. If long, rows are
            observations with three cols each: classification_id, tax_id, and count.

        Returns
        -------
        ClassificationsDataFrame
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

        if normalize is True or (normalize == 'auto' and rank is not None and self._field != 'abundance'):
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

        # additional data to copy into the ClassificationsDataFrame
        ocx_data = {
            'ocx_metadata': self.metadata.copy(),
            'ocx_rank': rank,
            'ocx_field': self._field,
            'ocx_taxonomy': self.taxonomy.copy(),
        }

        # generate long-format table
        if table_format == 'long':
            long_df = {
                'classification_id': [],
                'tax_id': [],
                self._field: []
            }

            for t_id in df:
                for c_id, count in df[t_id].iteritems():
                    long_df['classification_id'].append(c_id)
                    long_df['tax_id'].append(t_id)
                    long_df[self._field].append(count)

            results_df = ClassificationsDataFrame(long_df, **ocx_data)
        elif table_format == 'wide':
            results_df = ClassificationsDataFrame(df, **ocx_data)
        else:
            raise OneCodexException('table_format must be one of: long, wide')

        return results_df


@pd.api.extensions.register_dataframe_accessor('ocx')
class OneCodexAccessor(AnalysisMixin):
    def __init__(self, pandas_obj):
        # copy data from the ClassificationsDataFrame to a new instance of AnalysisMethods
        self.metadata = pandas_obj.ocx_metadata
        self.taxonomy = pandas_obj.ocx_taxonomy
        self._field = pandas_obj.ocx_field
        self._rank = pandas_obj.ocx_rank
        self._results = pandas_obj

        # prune back _taxonomy df to contain only taxa present in the ClassificationsDataFrame (and parents)
        tree = self.tree_build()
        tree = self.tree_prune_tax_ids(tree, self._results.keys())

        tax_ids_to_keep = [x.name for x in tree.traverse()]

        self.taxonomy = self.taxonomy.loc[tax_ids_to_keep]

        # similarly restrict _metadata df to contain only data relevant to samples currently in ClassificationsDataFrame
        self.metadata = self.metadata.loc[self._results.index]
