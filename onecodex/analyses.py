from __future__ import annotations

import warnings
from collections import Counter
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Optional, Union, Literal

from onecodex.exceptions import OneCodexException
from onecodex.lib.enums import (
    AbundanceMetric,
    AnalysisType,
    FunctionalAnnotations,
    FunctionalAnnotationsMetric,
    Metric,
    Rank,
)
from onecodex.stats import StatsMixin
from onecodex.viz import (
    VizBargraphMixin,
    VizDistanceMixin,
    VizFunctionalHeatmapMixin,
    VizHeatmapMixin,
    VizMetadataMixin,
    VizPCAMixin,
)

if TYPE_CHECKING:
    import pandas as pd


def _get_classification_ids_without_abundances(df):
    classification_ids_without_abundances = []
    for class_id, is_all_nan in df.isnull().all(1).items():
        if is_all_nan:
            classification_ids_without_abundances.append(class_id)
    return classification_ids_without_abundances


@dataclass
class MetadataFetchResults:
    # Transformed metadata
    df: pd.DataFrame
    # Map of original fields to renamed fields
    renamed_fields: dict[str | int | tuple[str, ...], str]
    # Set of *original* fields that matched taxonomy
    taxonomy_fields: set[str | int | tuple[str, ...]]


class AnalysisMixin(
    VizPCAMixin,
    VizHeatmapMixin,
    VizMetadataMixin,
    VizDistanceMixin,
    VizBargraphMixin,
    VizFunctionalHeatmapMixin,
    StatsMixin,
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
        from onecodex.dataframes import OneCodexAccessor

        if rank == Rank.Auto:
            # if we're an accessor for a ClassificationsDataFrame, use its _rank property
            if isinstance(self, OneCodexAccessor):
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

    def _metadata_fetch(
        self, metadata_fields, label=None, match_taxonomy=True, coerce_missing_composite_fields=True
    ) -> MetadataFetchResults:
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
        match_taxonomy : `bool`, optional
            Whether the returned metadata should resolve `metadata_fields` against taxonomy
            if they're not included in sample metadata. Defaults to true.
        coerce_missing_composite_fields : bool, optional
            If a field in `metadata_fields` is a tuple, by default any missing values will be
            coerced to strings when creating the composite field. For example, `None` will become
            `"None"` and `np.nan` will become `"nan"`. If `False`, the composite field will be
            `np.nan` if any of its elements contain missing data.

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
        Dataclass with the following attributes:

        `df`: `pandas.DataFrame`
            Columns are renamed (if applicable) metadata fields and rows are `Classifications.id`.
            Elements are transformed values. Not all metadata fields will have been renamed, but will
            be present in the below `dict` nonetheless.
        `renamed_fields`: `dict`
            Keys are metadata fields and values are renamed metadata fields. This can be used to map
            metadata fields which were passed to this function, to prettier names. For example, if
            'bacteroid' is passed, it will be matched with the Bacteroides genus and renamed to
            'Bacteroides (816)', which includes its taxon ID.
        `taxonomy_fields`: `set`
            Set of *original* fields that matched taxonomy (either by taxon ID or taxon name).

        """
        import numpy as np
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

        # keep track of which fields match taxonomy
        taxonomy_fields = set()

        for f in set([f for f in metadata_fields if f]):
            if isinstance(f, tuple):
                # joined categorical metadata
                for field in f:
                    if not isinstance(field, str):
                        raise OneCodexException(f"Metadata field name {field} must be of type str")

                    if field not in self.metadata:
                        raise OneCodexException(
                            f"Metadata field {field} not found. Choose from: {help_metadata}"
                        )

                    if not (
                        pd.api.types.is_bool_dtype(self.metadata[field])
                        or isinstance(self.metadata[field].dtype, pd.CategoricalDtype)
                        or pd.api.types.is_object_dtype(self.metadata[field])
                    ):
                        raise OneCodexException(
                            "When specifying multiple metadata fields, all must be categorical"
                        )

                # concatenate the columns together with underscores
                composite_field = "_".join(f)
                magic_metadata[composite_field] = ""

                if coerce_missing_composite_fields:
                    magic_metadata[composite_field] = (
                        magic_metadata[composite_field]
                        .str.cat([self.metadata[field].astype(str) for field in f], sep="_")
                        .str.lstrip("_")
                    )
                else:
                    # https://stackoverflow.com/a/47333556/3776794
                    magic_metadata[composite_field] = (
                        magic_metadata[composite_field]
                        .str.cat(
                            [
                                np.where(
                                    pd.isnull(self.metadata[field]),
                                    self.metadata[field],
                                    self.metadata[field].astype(str),
                                )
                                for field in f
                            ],
                            sep="_",
                            na_rep=None,
                        )
                        .str.lstrip("_")
                    )

                magic_fields[f] = composite_field
            else:
                str_f = str(f)

                if str_f == "Label":
                    magic_metadata[str_f] = self.metadata["filename"]
                    magic_fields[f] = str_f
                    if label is not None:
                        magic_metadata[str_f] = self._make_labels_by_item_id(self.metadata, label)
                elif str_f in self.metadata:
                    # exactly matches existing metadata field
                    magic_metadata[f] = self.metadata[str_f]
                    magic_fields[f] = str_f
                elif match_taxonomy and str_f in self._results.keys():
                    # is a tax_id
                    rank = self.taxonomy["rank"][str_f]
                    if Rank.has_value(rank):
                        tax_name = self.taxonomy["name"][str_f]

                        # report within-rank abundance
                        df = self.to_df(rank=rank)

                        renamed_field = "{} ({})".format(tax_name, str_f)
                        magic_metadata[renamed_field] = df[str_f]
                        magic_fields[f] = renamed_field
                        taxonomy_fields.add(f)
                    else:
                        # matched a non-canonical rank
                        magic_metadata[f] = None
                        magic_fields[f] = str_f
                elif match_taxonomy:
                    # try to match it up with a taxon name
                    hits = []

                    # don't both searching if the query is really short
                    if len(str_f) > 4:
                        for tax_id, rank, tax_name in zip(
                            self.taxonomy.index, self.taxonomy["rank"], self.taxonomy["name"]
                        ):
                            if not Rank.has_value(rank):
                                # don't match non-canonical ranks
                                continue

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
                        taxonomy_fields.add(f)
                    else:
                        # matched nothing
                        magic_metadata[f] = None
                        magic_fields[f] = str_f
                else:
                    # matched nothing
                    magic_metadata[f] = None
                    magic_fields[f] = str_f

        return MetadataFetchResults(
            df=magic_metadata, renamed_fields=magic_fields, taxonomy_fields=taxonomy_fields
        )

    @property
    def _classification_ids_without_abundances(self):
        """
        Provide list of classification ids for which there are no abundances calculated.

        This can be used in plotting functions that may want to exclude or represent samples or
        classifications in this category differently. Since plotting functions can be called on
        `SampleCollection`s or `OneCodexAccessor`s, storing this list is safest, before
        dataframe manipulation happens and information is possibly lost.
        """
        from onecodex.dataframes import OneCodexAccessor

        if not AbundanceMetric.has_value(self._metric):
            return []

        if isinstance(self, OneCodexAccessor):
            if self._ocx_classification_ids_without_abundances is None:
                # We rely on this list for accurate plot data, and it shouldn't ever be None, but if
                # it is None, we raise an exception to avoid generating misleading plots
                raise OneCodexException(
                    "Unable to fetch list of classification IDs without abundances"
                )
            return self._ocx_classification_ids_without_abundances
        return _get_classification_ids_without_abundances(self._results)

    def to_df(self, analysis_type=AnalysisType.Classification, **kwargs):
        """
        Transform Analyses of samples in a `SampleCollection` into tabular format.

        Parameters
        ----------
        analysis_type : {'classification', 'functional'}, optional
            The `analysis_type` to aggregate, corresponding to AnalysisJob.analysis_type
        kwargs : dict, optional
             Keyword arguments specific to the `analysis_type`; see each individual function definition

        .. seealso:: to_classification_df
        .. seealso:: to_functional_df
        """
        generate_df = {
            AnalysisType.Classification: self.to_classification_df,
            AnalysisType.Functional: self.to_functional_df,
        }
        return generate_df[AnalysisType(analysis_type)](**kwargs)

    def to_functional_df(
        self,
        annotation: FunctionalAnnotations = FunctionalAnnotations.Pathways,
        taxa_stratified: bool = True,
        metric: FunctionalAnnotationsMetric = FunctionalAnnotationsMetric.Coverage,
        fill_missing: bool = True,
        filler: Any = 0,
    ):
        """
        Generate a FunctionalDataFrame associated with functional analysis results.

        Parameters
        ----------
        annotation : {onecodex.lib.enum.FunctionalAnnotations, str}, optional
            Annotation data to return, defaults to `pathways`
        taxa_stratified : bool, optional
            Return taxonomically stratified data, defaults to `True`
        metric : {onecodex.lib.enum.FunctionalAnnotationsMetric, str}, optional
            Metric values to return
            {'coverage', 'abundance'} for annotation==FunctionalAnnotations.Pathways or
            {'rpk', 'cpm'} for other annotations, defaults to `coverage`
        fill_missing : bool, optional
            Fill np.nan values
        filler : float, optional
            Value with which to fill np.nans
        """
        from onecodex.dataframes import FunctionalDataFrame

        df, feature_name_map = self._functional_results(
            annotation=annotation,
            taxa_stratified=taxa_stratified,
            metric=metric,
            fill_missing=fill_missing,
            filler=filler,
        )
        return FunctionalDataFrame(
            df,
            ocx_metadata=self.metadata.copy(),
            ocx_functional_group=annotation,
            ocx_metric=metric,
            ocx_feature_name_map=feature_name_map,
        )

    def to_classification_df(
        self,
        rank: Rank = Rank.Auto,
        top_n: Optional[int] = None,
        threshold: Optional[float] = None,
        remove_zeros: bool = True,
        normalize: Union[Literal["auto"], bool] = "auto",
        table_format: Union[Literal["wide"] | Literal["long"]] = "wide",
        include_taxa_missing_rank: bool = False,
        fill_missing: bool = True,
        filler: Any = 0,
    ):
        """Generate a ClassificationsDataFrame, performing any specified transformations.

        Takes the ClassificationsDataFrame associated with these samples, or SampleCollection,
        does some filtering, and returns a ClassificationsDataFrame copy.

        Parameters
        ----------
        rank : {'auto', 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.
        top_n : `integer`, optional
            Return only the top N most abundant taxa.
        threshold : `float`, optional
            Return only taxa more abundant than this threshold in one or more samples.
        remove_zeros : `bool`, optional
            Do not return taxa that have zero abundance in every sample.
        normalize : {'auto', True, False}
            Convert read counts to relative abundances (each sample sums to 1.0). If data has
            already been normalized, passing ``normalize=False`` will raise an error. To generate
            denormalized data, please create a new SampleCollection with ``metric="readcount"`` or
            ``metric="readcount_w_children"``.
        table_format : {'long', 'wide'}
            If wide, rows are classifications, cols are taxa, elements are counts. If long, rows are
            observations with three cols each: classification_id, tax_id, and count.
        include_taxa_missing_rank : bool, optional
            Whether or not to include taxa that do not have a designated parent at `rank` (will be
            grouped into a "No <rank>" column).
        fill_missing : bool, optional
            Fill np.nan values
        filler : float, optional
            Value with which to fill np.nans

        Returns
        -------
        `ClassificationsDataFrame`
        """
        from onecodex.dataframes import ClassificationsDataFrame

        if include_taxa_missing_rank:
            if not rank:
                raise OneCodexException(
                    "`rank` must be specified when passing `include_taxa_missing_rank=True`."
                )

            if self._metric not in {Metric.ReadcountWChildren, Metric.AbundanceWChildren}:
                raise OneCodexException(
                    "`include_taxa_missing_rank` can only be used with `readcount_w_children` or "
                    "`abundance_w_children` metrics."
                )

        rank = self._get_auto_rank(rank)

        df = self._results.copy()
        if fill_missing:
            df = df.fillna(filler)

        # subset by taxa
        if rank:
            try:
                rank = Rank(rank)
            except ValueError:
                raise OneCodexException(f"Invalid rank: {rank}")

            if rank == Rank.Kingdom:
                warnings.warn(
                    "Did you mean to specify rank=kingdom? Use rank=superkingdom to see Bacteria, "
                    "Archaea and Eukaryota."
                )

            level = rank.level
            tax_ids_to_keep = []
            unclassified_tax_ids = set()

            for tax_id in df.keys():
                try:
                    tax_id_level = Rank(self.taxonomy["rank"][tax_id]).level
                except ValueError:
                    tax_id_level = None

                if tax_id_level == level:
                    tax_ids_to_keep.append(tax_id)
                elif (
                    include_taxa_missing_rank and tax_id_level is not None and tax_id_level < level
                ):
                    unclassified_tax_id = self._get_highest_unclassified_tax_id(tax_id, level)

                    if unclassified_tax_id is not None:
                        unclassified_tax_ids.add(unclassified_tax_id)

            if unclassified_tax_ids:
                no_level_name = f"No {rank.value}"
                df[no_level_name] = df[list(unclassified_tax_ids)].sum(axis=1)
                tax_ids_to_keep.append(no_level_name)

            if len(tax_ids_to_keep) == 0:
                raise OneCodexException(f"No taxa kept--is rank ({rank.value}) correct?")

            df = df.loc[:, tax_ids_to_keep]

        # normalize
        if normalize is False and self._guess_normalized():
            raise OneCodexException(
                f"Data has already been normalized. To generate denormalized data, please create a "
                f"new SampleCollection with `metric={Metric.Readcount.value!r}` or "
                f"`metric={Metric.ReadcountWChildren.value!r}`."
            )

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
            "ocx_classification_ids_without_abundances": self._classification_ids_without_abundances,
        }

        # generate long-format table
        if table_format == "long":
            pretty_metric_name = self._make_pretty_metric_name(self._metric, normalize)

            long_df = {"classification_id": [], "tax_id": [], pretty_metric_name: []}

            for t_id in df:
                for c_id, count in df[t_id].items():
                    long_df["classification_id"].append(c_id)
                    long_df["tax_id"].append(t_id)
                    long_df[pretty_metric_name].append(count)

            results_df = ClassificationsDataFrame(long_df, **ocx_data)
        elif table_format == "wide":
            results_df = ClassificationsDataFrame(df, **ocx_data)
        else:
            raise OneCodexException("table_format must be one of: long, wide")

        return results_df

    def _get_highest_unclassified_tax_id(self, tax_id, level):
        # try to determine if child nodes have parent nodes at the specified level. we use this to
        # make counts of branches that "disappear" at different levels so we can insert an extra
        # abundance for them (e.g. "No genus")
        curr_tax_id = tax_id
        highest_tax_id = tax_id

        while curr_tax_id is not None:
            try:
                curr_level = Rank(self.taxonomy["rank"][curr_tax_id]).level
            except ValueError:
                curr_level = None

            if curr_level is not None:
                if curr_level > level:
                    break

                if curr_level == level:
                    # there's something at the level we're interested in, so this is actually
                    # classified and we don't have to keep track of it
                    return None

                # the highest tax id with a level definitely under the level we're classifying to
                highest_tax_id = curr_tax_id

            curr_tax_id = self.taxonomy["parent_tax_id"][curr_tax_id]

        return highest_tax_id

    @staticmethod
    def _make_labels_by_item_id(metadata, label):
        """Make/Extract labels from metadata pandas dataframe.

        Parameters
        ----------
        metadata : `pandas.DataFrame`
        label : `str` or `callable`

        Returns
        -------
        `dict`
            Keys are from metadata.index. Values are generated labels.
        """
        import pandas as pd

        raw_result = {}

        if isinstance(label, str):
            if label in metadata.columns:
                raw_result = dict(metadata[label].astype(str).items())
            else:
                raise OneCodexException(
                    "Label field {} not found. Choose from: {}".format(
                        label, ", ".join(metadata.keys())
                    )
                )
        elif callable(label):
            for item_id, item_meta in metadata.to_dict(orient="index").items():
                item_label = label(item_meta)
                if not isinstance(item_label, str):
                    wrong_type = type(item_label).__name__
                    raise OneCodexException(
                        "Expected string from label function, got: {}".format(wrong_type)
                    )
                raw_result[item_id] = item_label

        elif label is not None:
            wrong_type = type(label).__name__
            raise OneCodexException(
                "Expected string or callable for label, got: {}".format(wrong_type)
            )

        # add an incremented number to duplicate labels (e.g., same filename)
        duplicates_counter = Counter(raw_result.values())
        duplicated_labels = [label for label, n in duplicates_counter.items() if n > 1]
        indexing = {x: 1 for x in duplicated_labels}
        result = {}
        for item_id, label in raw_result.items():
            if label in indexing:
                result[item_id] = "{} ({})".format(label, indexing[label])
                indexing[label] += 1
            else:
                result[item_id] = label

        return pd.Series(result)

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
