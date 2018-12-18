import pandas as pd
import warnings

import six

from onecodex.exceptions import OneCodexException
from onecodex.models import Analyses, Classifications, Samples, Metadata
from onecodex.viz import VizPCAMixin, VizHeatmapMixin, VizMetadataMixin, VizDistanceMixin
from onecodex.distance import DistanceMixin


# force persistence of our additional taxonomy and metadata dataframe properties
class ResultsDataFrame(pd.DataFrame):
    _metadata = ['ocx_field', 'ocx_taxonomy', 'ocx_metadata']

    def __init__(self, data=None, index=None, columns=None, dtype=None, copy=False, ocx_data={}):
        try:
            self.ocx_field = ocx_data['ocx_field']
            self.ocx_taxonomy = ocx_data['ocx_taxonomy']
            self.ocx_metadata = ocx_data['ocx_metadata']
        except KeyError:
            pass

        pd.DataFrame.__init__(self, data=data, index=index, columns=columns, dtype=dtype, copy=copy)

    @property
    def _constructor(self):
        from functools import partial

        ocx_data = {
            'ocx_field': self.ocx_field,
            'ocx_taxonomy': self.ocx_taxonomy,
            'ocx_metadata': self.ocx_metadata
        }

        return partial(ResultsDataFrame, ocx_data=ocx_data)

    @property
    def _constructor_sliced(self):
        return ResultsSeries


class ResultsSeries(pd.Series):
    _metadata = ['ocx_field', 'ocx_taxonomy', 'ocx_metadata']

    @property
    def _constructor(self):
        return ResultsSeries

    @property
    def _constructor_expanddim(self):
        return ResultsDataFrame


class SampleCollection(VizPCAMixin, VizHeatmapMixin, VizMetadataMixin, DistanceMixin, VizDistanceMixin):
    _immutable = ['_immutable', 'field', '_analyses', '_results', '_metadata', '_taxonomy', '_classifications']

    def __init__(self, analyses, **kwargs):
        self.__setattr__('_analyses', analyses, force=True)

        self.magic_classification_fetch(**kwargs)
        self.collate_metadata(**kwargs)
        self.collate_results(**kwargs)

    def __setattr__(self, name, value, force=False):
        if name in self._immutable and not force:
            raise OneCodexException('{} is an immutable property and can not be changed'.format(name))

        super(SampleCollection, self).__setattr__(name, value)

    def _get_auto_rank(self, rank):
        if rank == 'auto':
            if self.field == 'abundance':
                return 'species'
            else:
                return 'genus'
        else:
            return rank

    def magic_classification_fetch(self, skip_missing=True, **kwargs):
        """Turns a list of objects associated with a classification results into a list of
        Classifications objects.

        analyses (list) -- list of Samples, Classifications, or Analyses objects
        skip_missing (bool) -- if an analysis was not successful, exclude it and keep going
        warn (bool) -- issue warnings
        """

        fetched = []

        for a in self._analyses:
            if isinstance(a, Samples):
                c = a.primary_classification
            elif isinstance(a, Classifications):
                c = a
            elif isinstance(a, Analyses):
                if a.analysis_type != 'classification':
                    raise OneCodexException('{} is not a classification'.format(a.id))
                c = Classifications(a._resource._client.Classifications.fetch(a.id))

            if skip_missing and not c.success:
                warnings.warn('Classification {} not successful. Skipping.'.format(c.id))
                continue

            fetched.append(c)

        self.__setattr__('_classifications', fetched, force=True)

    def collate_metadata(self, label=None, **kwargs):
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

    def collate_results(self, field='auto', **kwargs):
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

        # # if field is 'auto' we will try to use abundances but fall back on readcount_w_children
        # if field == 'auto':
        #     c_with_abundance = []

        #     # does every classification result given have 'abundance' field?
        #     for c_idx, c in enumerate(self._classifications):
        #         result = c.results()['table']

        #         for d in result[:10]:
        #             if 'abundance' not in d:
        #                 c_with_abundance.append(False)
        #                 break
        #         else:
        #             c_with_abundance.append(True)

        #     if all(c_with_abundance):
        #         field = 'abundance'
        #     else:
        #         if any(c_with_abundance):
        #             warnings.warn('Be aware, you are mixing shotgun and 16S classification results!')

        #         field = 'readcount_w_children'

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
                     .set_index('tax_id') \

        self.__setattr__('_results', df, force=True)
        self.__setattr__('_taxonomy', tax_info, force=True)

    def results(self, rank=None, top_n=None, threshold=None,
                remove_zeros=True, normalize='auto',
                table_format='wide', **kwargs):
        """
        Filtering options (performed consecutively in the order given below):
            rank (None | 'kingdom' | 'phylum' | 'class' | 'order' | 'family' | 'genus' | 'species')
                - None: include all ranks
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

        results_df.ocx_field = self.field
        results_df.ocx_taxonomy = self.taxonomy()
        results_df.ocx_metadata = self.metadata()

        return results_df

    def metadata(self):
        return self._metadata.copy()

    def taxonomy(self):
        return self._taxonomy.copy()

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

                composite_field = '_'.join(f)
                magic_metadata[composite_field] = ''
                magic_metadata[composite_field] = magic_metadata[composite_field].str.cat(
                    [self._metadata[field].astype(str) for field in f], sep='_'
                ).str.lstrip('_')
                magic_fields[f] = composite_field
            else:
                f = str(f)

                if f == 'Label':
                    # it's our magic keyword that means _display_name
                    magic_metadata[f] = self._metadata['_display_name']
                    magic_fields[f] = f
                elif f in self._metadata:
                    # exactly matches existing metadata field
                    magic_metadata[f] = self._metadata[f]
                    magic_fields[f] = f
                elif f in self._results.keys():
                    # is a tax_id
                    tax_name = self._taxonomy['name'][f]

                    # report within-rank abundance
                    df = self.results(rank=self._taxonomy['rank'][f])

                    renamed_field = '{} ({})'.format(tax_name, f)
                    magic_metadata[renamed_field] = df[f]
                    magic_fields[f] = renamed_field
                else:
                    # try to match it up with a taxon name
                    hits = []

                    for tax_id, tax_name in zip(self._taxonomy.index, self._taxonomy['name']):
                        # if it's an exact match, use that and skip the rest
                        if f.lower() == tax_name.lower():
                            hits = [(tax_id, tax_name)]
                            break
                        # otherwise, keep trying to match
                        elif f.lower() in tax_name.lower():
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
                            'Field or taxon {} not found. Choose from: {}'.format(f, help_metadata)
                        )

        return magic_metadata, magic_fields


@pd.api.extensions.register_dataframe_accessor('ocx')
class OneCodexAccessor(SampleCollection):
    def __init__(self, pandas_obj):
        super(OneCodexAccessor, self).__setattr__('_results', pandas_obj, force=True)
        super(OneCodexAccessor, self).__setattr__('field', pandas_obj.ocx_field, force=True)
        super(OneCodexAccessor, self).__setattr__('_taxonomy', pandas_obj.ocx_taxonomy, force=True)
        super(OneCodexAccessor, self).__setattr__('_metadata', pandas_obj.ocx_metadata, force=True)
