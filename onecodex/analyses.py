import six
import warnings

from onecodex.exceptions import OneCodexException
from onecodex.lib.enums import AbundanceMetric, Rank
from onecodex.viz import (
    VizPCAMixin,
    VizHeatmapMixin,
    VizMetadataMixin,
    VizDistanceMixin,
    VizBargraphMixin,
)


class AnalysisMixin(
    VizPCAMixin, VizHeatmapMixin, VizMetadataMixin, VizDistanceMixin, VizBargraphMixin
):
    """Contains methods for analyzing Classifications results.

    Notes
    -----
    Three DataFrames are required by most methods: collated counts, collated metadata, and taxonomy.
    This data is obtained from either a `ClassificationsDataFrame` or a `SampleCollection`. Both
    classes use this mixin. `AnalysisMixin` pulls additional methods in from `onecodex.distance`,
    `onecodex.taxonomy`, and `onecodex.viz`.
    """

    def _get_auto_rank(self, rank):
        """Attempt to figure out what rank we should use for analyses."""

        if rank == Rank.Auto:
            # if we're an accessor for a ClassificationsDataFrame, use its _rank property
            if self.__class__.__name__ == "OneCodexAccessor":
                return self._rank

            if AbundanceMetric.has_value(self._metric) or self._is_metagenomic:
                return Rank.Species
            else:
                return Rank.Genus
        else:
            return rank

    def _guess_normalized(self):
        """Return True if the collated counts in `self._results` appear to be normalized.

        Notes
        -----
        It's possible that the _results df has already been normalized, which can cause some
        methods to fail. This method lets us guess whether that's true and act accordingly.
        """
        return (
            getattr(self, "_normalized", False)
            or AbundanceMetric.has_value(self._metric)
            or bool((self._results.sum(axis=1).round(4) == 1.0).all())
        )  # noqa

    def _metadata_fetch(self, metadata_fields, label=None):
        """Fetch and transform given metadata fields from `self.metadata`.

        Takes a list of metadata fields, some of which can contain taxon names or taxon IDs, and
        returns a DataFrame with transformed data that can be used for plotting.

        Parameters
        ----------
        metadata_fields : `list` of `string`
            A list of metadata fields, taxon names, or taxon IDs to fetch and transform for display.
        label : `string` or `callable`, optional
            A metadata field (or function) used to label each analysis. If passing a function, a
            dict containing the metadata for each analysis is passed as the first and only
            positional argument. The callable function must return a string.

            If this argument is not given, and "Label" is in `metadata_fields`, "Label" will be set
            to the filename associated with an analysis.

        Notes
        -----
        Taxon names and IDs are transformed into the relative abundances of those taxa within their
        own rank. For example, 'Bacteroides' will return the relative abundances of 'Bacteroides'
        among all taxa of rank genus. Taxon IDs are stored as strings in `ClassificationsDataFrame`
        and are coerced to strings if integers are given.

        Metadata fields are returned as is, from the `self.metadata` DataFrame. If multiple metadata
        fields are specified in a tuple, their values are joined as strings separated by underscore.
        Multiple metadata fields in a tuple must both be categorical. That is, a numerical field and
        boolean can not be joined, or the result would be something like '87.4_True'.

        Returns
        -------
        `pandas.DataFrame`
            Columns are renamed (if applicable) metadata fields and rows are `Classifications.id`.
            Elements are transformed values. Not all metadata fields will have been renamed, but will
            be present in the below `dict` nonetheless.
        `dict`
            Keys are metadata fields and values are renamed metadata fields. This can be used to map
            metadata fields which were passed to this function, to prettier names. For example, if
            'bacteroid' is passed, it will be matched with the Bacteroides genus and renamed to
            'Bacteroides (816)', which includes its taxon ID.
        """
        import pandas as pd

        help_metadata = ", ".join(self.metadata.keys())
        magic_metadata = pd.DataFrame({"classification_id": self._results.index}).set_index(
            "classification_id"
        )

        # if user passed label kwarg but didn't put "Label" in the fields, assume the user wants
        # that field added
        if label is not None and "Label" not in metadata_fields:
            metadata_fields.append("Label")

        # if we magically rename fields, keep track
        magic_fields = {}

        for f in set([f for f in metadata_fields if f]):
            if isinstance(f, tuple):
                # joined categorical metadata
                for field in f:
                    if field not in self.metadata:
                        raise OneCodexException(
                            "Metric {} not found. Choose from: {}".format(field, help_metadata)
                        )

                    if not (
                        pd.api.types.is_bool_dtype(self.metadata[field])
                        or pd.api.types.is_categorical_dtype(self.metadata[field])  # noqa
                        or pd.api.types.is_object_dtype(self.metadata[field])  # noqa
                    ):
                        raise OneCodexException(
                            "When specifying multiple metadata fields, all must be categorical"
                        )

                # concatenate the columns together with underscores
                composite_field = "_".join(f)
                magic_metadata[composite_field] = ""
                magic_metadata[composite_field] = (
                    magic_metadata[composite_field]
                    .str.cat([self.metadata[field].astype(str) for field in f], sep="_")
                    .str.lstrip("_")
                )
                magic_fields[f] = composite_field
            else:
                str_f = str(f)

                if str_f == "Label":
                    magic_metadata[str_f] = self.metadata["filename"]
                    magic_fields[f] = str_f

                    if isinstance(label, six.string_types):
                        if label in self.metadata.columns:
                            magic_metadata[str_f] = self.metadata[label].astype(str)
                        else:
                            raise OneCodexException(
                                "Label field {} not found. Choose from: {}".format(
                                    label, help_metadata
                                )
                            )
                    elif callable(label):
                        for classification_id, metadata in self.metadata.to_dict(
                            orient="index"
                        ).items():
                            c_id_label = label(metadata)

                            if not isinstance(c_id_label, six.string_types):
                                raise OneCodexException(
                                    "Expected string from label function, got: {}".format(
                                        type(c_id_label).__name__
                                    )
                                )

                            magic_metadata.loc[classification_id, "Label"] = c_id_label
                    elif label is not None:
                        raise OneCodexException(
                            "Expected string or callable for label, got: {}".format(
                                type(label).__name__
                            )
                        )

                    # add an incremented number to duplicate labels (e.g., same filename)
                    duplicate_labels = (
                        magic_metadata[str_f]
                        .where(magic_metadata[str_f].duplicated(keep=False))
                        .dropna()
                    )

                    if not duplicate_labels.empty:
                        duplicate_counts = {label: 1 for label in duplicate_labels}

                        for c_id in magic_metadata.index:
                            label = magic_metadata[str_f][c_id]

                            if duplicate_labels.isin([label]).any():
                                magic_metadata[str_f][c_id] = "{} ({})".format(
                                    label, duplicate_counts[label]
                                )
                                duplicate_counts[label] += 1
                elif str_f in self.metadata:
                    # exactly matches existing metadata field
                    magic_metadata[f] = self.metadata[str_f]
                    magic_fields[f] = str_f
                elif str_f in self._results.keys():
                    # is a tax_id
                    tax_name = self.taxonomy["name"][str_f]

                    # report within-rank abundance
                    df = self.to_df(rank=self.taxonomy["rank"][str_f])

                    renamed_field = "{} ({})".format(tax_name, str_f)
                    magic_metadata[renamed_field] = df[str_f]
                    magic_fields[f] = renamed_field
                else:
                    # try to match it up with a taxon name
                    hits = []

                    # don't both searching if the query is really short
                    if len(str_f) > 4:
                        for tax_id, tax_name in zip(self.taxonomy.index, self.taxonomy["name"]):
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
                        df = self.to_df(rank=self.taxonomy["rank"][hits[0][0]])

                        renamed_field = "{} ({})".format(hits[0][1], hits[0][0])
                        magic_metadata[renamed_field] = df[hits[0][0]]
                        magic_fields[f] = renamed_field
                    else:
                        # matched nothing
                        raise OneCodexException(
                            "Metric or taxon {} not found. Choose from: {}".format(
                                str_f, help_metadata
                            )
                        )

        return magic_metadata, magic_fields

    def to_df(
        self,
        rank=Rank.Auto,
        top_n=None,
        threshold=None,
        remove_zeros=True,
        normalize="auto",
        table_format="wide",
    ):
        """Generate a ClassificationDataFrame, performing any specified transformations.

        Takes the ClassificationsDataFrame associated with these samples, or SampleCollection,
        does some filtering, and returns a ClassificationsDataFrame copy.

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
        `ClassificationsDataFrame`
        """
        from onecodex.dataframes import ClassificationsDataFrame

        rank = self._get_auto_rank(rank)
        df = self._results.copy()

        # subset by taxa
        if rank:
            if rank == "kingdom":
                warnings.warn(
                    "Did you mean to specify rank=kingdom? Use rank=superkingdom to see Bacteria, "
                    "Archaea and Eukaryota."
                )

            tax_ids_to_keep = []

            for tax_id in df.keys():
                if self.taxonomy["rank"][tax_id] == rank:
                    tax_ids_to_keep.append(tax_id)

            if len(tax_ids_to_keep) == 0:
                raise OneCodexException("No taxa kept--is rank ({}) correct?".format(rank))

            df = df.loc[:, tax_ids_to_keep]

        # normalize
        if normalize is False and self._guess_normalized():
            raise OneCodexException("Data has already been normalized and this can not be undone.")

        if normalize is True or (normalize == "auto" and rank):
            if not self._guess_normalized():
                # Replace nans with zeros for samples that have a total abundance of zero.
                df = df.div(df.sum(axis=1), axis=0).fillna(0.0)

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
            "ocx_metadata": self.metadata.copy(),
            "ocx_rank": rank,
            "ocx_metric": self._metric,
            "ocx_taxonomy": self.taxonomy.copy(),
            "ocx_normalized": normalize,
        }

        # generate long-format table
        if table_format == "long":
            pretty_metric_name = self._make_pretty_metric_name(self._metric, normalize)

            long_df = {"classification_id": [], "tax_id": [], pretty_metric_name: []}

            for t_id in df:
                for c_id, count in df[t_id].iteritems():
                    long_df["classification_id"].append(c_id)
                    long_df["tax_id"].append(t_id)
                    long_df[pretty_metric_name].append(count)

            results_df = ClassificationsDataFrame(long_df, **ocx_data)
        elif table_format == "wide":
            results_df = ClassificationsDataFrame(df, **ocx_data)
        else:
            raise OneCodexException("table_format must be one of: long, wide")

        return results_df

    @staticmethod
    def _make_pretty_metric_name(metric, normalized):
        if AbundanceMetric.has_value(metric):
            return "Relative Abundance"
        if normalized:
            return "Reads (Normalized)"
        else:
            return "Reads"

    @property
    def metric(self):
        return self._make_pretty_metric_name(self._metric, self._guess_normalized())
