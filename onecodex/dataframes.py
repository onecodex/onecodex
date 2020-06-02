import pandas as pd

from onecodex.analyses import AnalysisMixin


class ClassificationsDataFrame(pd.DataFrame):
    """A DataFrame containing additional One Codex metadata.

    A subclassed `pandas.DataFrame` containing additional metadata pertinent to analysis of
    One Codex Classifications results. These fields, once part of the DataFrame, will no longer be
    updated when the contents of the associated `SampleCollection` change. In comparison, the
    corresponding attributes `_rank`, `_metric`, `taxonomy` and `metadata` in a `SampleCollection`
    are re-generated whenever members of the `SampleCollection` are added or removed.

    Methods from `AnalysisMixin`, such as `to_df`, are available via the `ocx` namespace. For
    example, `ClassificationsDataFrame().ocx.to_df()`.

    Parameters
    ----------
    ocx_rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
        Analysis was restricted to abundances of taxa at the specified level.

    ocx_metric : {'readcount_w_children', 'readcount', 'abundance'}
        Which metric was used for the abundance/count of a particular taxon in a sample.

        - 'readcount_w_children': total reads of this taxon and all its descendants
        - 'readcount': total reads of this taxon
        - 'abundance': genome size-normalized relative abundances, from shotgun sequencing

    ocx_metadata : `pandas.DataFrame`
        A DataFrame containing collated metadata fields for all samples in this analysis.

    ocx_taxonomy : `pandas.DataFrame`
        A DataFrame containing taxonomy information (i.e., id, name, rank, parent) for all taxa
        referenced in this analysis.

    ocx_normalized : `bool`
        Whether the results in this DataFrame were normalized, each sample summing to 1.0.
    """

    _metadata = ["ocx_rank", "ocx_metric", "ocx_taxonomy", "ocx_metadata", "ocx_normalized"]

    def __init__(
        self,
        data=None,
        index=None,
        columns=None,
        dtype=None,
        copy=False,
        ocx_rank=None,
        ocx_metric=None,
        ocx_taxonomy=None,
        ocx_metadata=None,
        ocx_normalized=None,
    ):
        self.ocx_rank = ocx_rank
        self.ocx_metric = ocx_metric
        self.ocx_taxonomy = ocx_taxonomy
        self.ocx_metadata = ocx_metadata
        self.ocx_normalized = ocx_normalized

        pd.DataFrame.__init__(self, data=data, index=index, columns=columns, dtype=dtype, copy=copy)

    @property
    def _constructor(self):
        return ClassificationsDataFrame

    @property
    def _constructor_sliced(self):
        return ClassificationsSeries

    def to_html(self, *args, **kwargs):
        classes = kwargs.pop("classes", [])
        if isinstance(classes, str):
            classes = [classes]
        classes.append("ocx_classifications_df")
        kwargs["classes"] = classes
        kwargs["float_format"] = "%0.3f"
        kwargs["max_rows"] = 15
        kwargs["max_cols"] = 10

        # round abundances to avoid long trails of zeros, and sort taxa in order of abundance
        if "classification_id" in self.columns and "tax_id" in self.columns:
            # long format
            field_name = set(self.columns).difference({"classification_id", "tax_id"}).pop()
            df = self.copy()
            df[field_name] = df[field_name].round(6)
            df = df.sort_values(field_name, ascending=False)
        else:
            # wide format
            df = self.round(6).reindex(columns=self.sum().sort_values(ascending=False).index)

        return super(ClassificationsDataFrame, df).to_html(*args, **kwargs)


class ClassificationsSeries(pd.Series):
    """A Series containing additional One Codex metadata.

    A subclassed `pandas.Series` containing additional metadata pertinent to analysis of
    One Codex Classifications results. See the docstring for `ClassificationsDataFrame`.
    """

    # 'name' is a piece of metadata specified by pd.Series--it's not ours
    _metadata = ["name", "ocx_rank", "ocx_metric", "ocx_taxonomy", "ocx_metadata", "ocx_normalized"]

    def __init__(
        self,
        data=None,
        index=None,
        dtype=None,
        name=None,
        copy=False,
        fastpath=False,
        ocx_rank=None,
        ocx_metric=None,
        ocx_taxonomy=None,
        ocx_metadata=None,
        ocx_normalized=None,
    ):
        self.ocx_rank = ocx_rank
        self.ocx_metric = ocx_metric
        self.ocx_taxonomy = ocx_taxonomy
        self.ocx_metadata = ocx_metadata
        self.ocx_normalized = ocx_normalized

        pd.Series.__init__(
            self, data=data, index=index, dtype=dtype, name=name, copy=copy, fastpath=fastpath
        )

    @property
    def _constructor(self):
        return ClassificationsSeries

    @property
    def _constructor_expanddim(self):
        return ClassificationsDataFrame


@pd.api.extensions.register_dataframe_accessor("ocx")
class OneCodexAccessor(AnalysisMixin):
    """A pandas accessor object with One Codex methods and metadata.

    Accessor object alllowing access of `AnalysisMixin` methods from the 'ocx' namespace of a
    `ClassificationsDataFrame`.

    Notes
    -----
    When instantiated, the accessor will prune the taxonomic tree back to contain only taxa
    referenced in the classification results (i.e., self._results). Similarly, metadata is sliced
    such that it contains only those `Classifications.id` in the results. This is because users may
    filter or modify the classification results to remove classification results (i.e., rows) or
    taxa (i.e., cols) from the `ClassificationsDataFrame` before accessing this namespace.
    """

    def __init__(self, pandas_obj):
        # copy data from the ClassificationsDataFrame to a new instance of AnalysisMethods
        self.metadata = pandas_obj.ocx_metadata
        self.taxonomy = pandas_obj.ocx_taxonomy
        self._normalized = pandas_obj.ocx_normalized
        self._metric = pandas_obj.ocx_metric
        self._rank = pandas_obj.ocx_rank
        self._results = pandas_obj

        # prune back _taxonomy df to contain only taxa and parents in the ClassificationsDataFrame
        tree = self.tree_build()

        if "classification_id" in self._results.columns and "tax_id" in self._results.columns:
            tree = self.tree_prune_tax_ids(tree, self._results["tax_id"].drop_duplicates())

            # similarly restrict _metadata df to contain only data relevant to samples in the df
            self.metadata = self.metadata.loc[self._results["classification_id"].drop_duplicates()]
        else:
            tree = self.tree_prune_tax_ids(tree, self._results.columns)
            self.metadata = self.metadata.loc[self._results.index]

        tax_ids_to_keep = [x.name for x in tree.traverse()]
        self.taxonomy = self.taxonomy.loc[tax_ids_to_keep]
