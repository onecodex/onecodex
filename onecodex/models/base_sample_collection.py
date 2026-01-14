from __future__ import annotations

import json
import warnings
from collections import Counter, OrderedDict, defaultdict
from collections.abc import MutableSequence
from dataclasses import dataclass
from datetime import datetime
from functools import cached_property, lru_cache
from typing import TYPE_CHECKING, Any, Literal, Optional

from onecodex.exceptions import OneCodexException
from onecodex.lib.enums import (
    AbundanceMetric,
    AnalysisType,
    FunctionalAnnotations,
    FunctionalAnnotationsMetric,
    Metric,
    Rank,
)
from onecodex.models.analysis import Classifications, FunctionalProfiles
from onecodex.models.sample import Samples

if TYPE_CHECKING:
    import pandas as pd

    from onecodex.dataframes import ClassificationsDataFrame
    from onecodex.models import Classifications, Jobs, Samples


@dataclass
class MetadataFetchResults:
    # Transformed metadata
    df: pd.DataFrame
    # Map of original fields to renamed fields
    renamed_fields: dict[str | int | tuple[str, ...], str]
    # Set of *original* fields that matched taxonomy
    taxonomy_fields: set[str | int | tuple[str, ...]]


CANONICAL_RANKS = (
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)


class BaseSampleCollection(
    MutableSequence,
):
    """
    Contains all of the methods used to fetch and normalize sample and analysis results data.

    BaseSampleCollection exists so that the source of methods used by the VizMixins can be resolved.
    The VizMixins are combined by inheritance in onecodex.models.SampleCollection.
    """

    def __hash__(self):
        """To enable pickling (for lru_cache)."""
        return hash(tuple(r.id for r in self._res_list))

    def _parse_classification_config_args(self, metric: Metric, rank: Rank) -> tuple[Metric, Rank]:
        """Validate arguments and determine defaults where necessary."""

        if rank is None:
            raise OneCodexException("Please specify a rank or 'auto' to choose automatically")

        if metric in (None, Metric.Auto):
            metric = self.automatic_metric

        if rank == Rank.Auto:
            rank = self.automatic_rank(metric=metric)

        return Metric.from_value(metric), Rank.from_value(rank)

    def __init__(
        self,
        objects: list[Samples] | list[Classifications],
        job: Jobs | None = None,
        skip_missing: bool = True,
    ):
        """Instantiate a new SampleCollection containing `Samples` or `Classifications` objects.

        Parameters
        ----------
        objects : list
            A list of `onecodex.models.Samples` or `onecodex.models.Classifications` objects
            which will be processed into a `SampleCollection`

        skip_missing : bool, optional
            If an analysis was not successful, exclude it, warn, and keep going

        Examples
        --------
        Given a list of Samples, create a new SampleCollection using abundances:

            samples = [sample1, sample2, sample3]
            collection = SampleCollection(samples)

        Notes
        -----
        To provide access to the list-like API of `ResourceList`, must also accept a list of
        unwrapped potion resources and a One Codex model.
        """

        if not all(
            [isinstance(obj, Samples) or isinstance(obj, Classifications) for obj in objects]
        ):
            raise OneCodexException(
                "SampleCollection can only contain One Codex Samples or Classifications objects"
            )

        # are they all the same model?
        object_classes = [type(obj) for obj in objects]

        if len(set(object_classes)) > 1:
            raise OneCodexException(
                "SampleCollection can contain Samples or Classifications, but not both"
            )

        model = objects[0].__class__ if len(objects) > 0 else Samples

        self._kwargs = {
            "skip_missing": skip_missing,
            "job": job,
        }
        self._oc_model = model
        self._res_list = objects

    def _check_valid_resource(self, other, check_for_dupes=True):
        if not isinstance(other, list):
            other = [other]

        other_ids = []
        for o in other:
            if not isinstance(o, self._oc_model):
                raise ValueError(
                    "Expected object of type '{}', got '{}'".format(
                        self._oc_model.__name__, type(o).__name__
                    )
                )

            other_ids.append(o.id)

        if check_for_dupes:
            # duplicates are not allowed
            self_ids = [s.id for s in self._res_list]

            if len(set(self_ids + other_ids)) != len(self_ids + other_ids):
                raise OneCodexException(
                    "{} cannot contain duplicate objects".format(self.__class__.__name__)
                )

    def __eq__(self, other):
        return all(id(x) == id(y) for x, y in zip(self, other))

    def __contains__(self, other):
        return self._res_list.__contains__(other)

    @property
    def __repr__(self):
        return self._res_list.__repr__

    @property
    def __len__(self):
        return self._res_list.__len__

    def __getitem__(self, x):
        wrapped = self._res_list[x]
        if isinstance(wrapped, list):
            return self.__class__(self._res_list[x], **self._kwargs)
        else:
            return wrapped

    def __setitem__(self, k, v):
        self._check_valid_resource(v)
        self._res_list[k] = v

    def __delitem__(self, x):
        del self._res_list[x]

    @property
    def __iter__(self):
        return self._res_list.__iter__

    @property
    def __reversed__(self):
        return self._res_list.__reversed__

    def __add__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError(
                'can only concatenate {} (not "{}") to {}'.format(
                    self.__class__.__name__, type(other), self.__class__.__name__
                )
            )
        new_obj = self.copy()
        new_obj.extend(other._res_list)
        return new_obj

    def append(self, x):
        self._check_valid_resource(x)
        self._res_list.append(x)

    def clear(self):
        self._res_list.clear()

    def copy(self):
        new_obj = self.__class__(self._res_list[:], self._oc_model, self._kwargs)
        return new_obj

    def count(self, x):
        # assume that ResourceList objects are identical if they share the same underlying resource
        self._check_valid_resource(x, check_for_dupes=False)
        n = 0
        for res_obj in self._res_list:
            if res_obj == x:
                n += 1
        return n

    def extend(self, iterable):
        self._check_valid_resource(iterable)
        self._res_list.extend([x for x in iterable])

    def index(self, x):
        # assume that ResourceList objects are identical if they share the same underlying resource
        self._check_valid_resource(x, check_for_dupes=False)
        for res_obj_idx, res_obj in enumerate(self._res_list):
            if res_obj == x:
                return res_obj_idx
        raise ValueError("{} is not in list".format(x))

    def insert(self, idx, x):
        self._check_valid_resource(x)
        self._res_list.insert(idx, x)

    def pop(self):
        return self._res_list.pop()

    def remove(self, x):
        del self._res_list[self.index(x)]

    def filter(self, filter_func):
        """Return a new `SampleCollection` containing only samples meeting the filter criteria.

        Will pass any kwargs (e.g., `metric` or `skip_missing`) used when instantiating the current class
        on to the new `SampleCollection` that is returned.

        Parameters
        ----------
        filter_func : callable
            A function that will be evaluated on every object in the collection. The function must
            return a `bool`. If True, the object will be kept. If False, it will be removed from the
            SampleCollection that is returned.

        Returns
        -------
        `onecodex.models.SampleCollection` containing only objects `filter_func` returned True on.

        Examples
        --------
        Generate a new collection of Samples that have a specific filename extension:

        >>> new_collection = samples.filter(lambda s: s.filename.endswith('.fastq.gz'))
        """
        if callable(filter_func):
            return self.__class__([obj for obj in self if filter_func(obj)], **self._kwargs)
        else:
            raise OneCodexException(
                "Please pass a function to filter: {}".format(type(filter_func).__name__)
            )

    @property
    def _classifications_from_res_list(self) -> list[Classifications | None]:
        from onecodex.models import Classifications, Samples

        classifications = []
        for obj in self._res_list:
            if isinstance(obj, Samples):
                classification = obj.primary_classification
            elif isinstance(obj, Classifications):
                classification = obj
            else:
                raise OneCodexException(
                    f"Objects in SampleCollection must be one of: Classifications, Samples, got {obj} {type(obj)}"
                )

            if not classification:
                msg = f"Classification not found for sample {obj.id}."
                if self._kwargs["skip_missing"]:
                    warnings.warn(msg + " Skipping.")
                    continue
                else:
                    raise OneCodexException(msg)

            classifications.append(classification)

        return classifications

    @property
    def _classifications(self) -> list[Classifications]:
        """Transform a list of Samples or Classifications into a list of Classifications objects.

        Parameters
        ----------
        skip_missing : bool
            If an analysis was not successful, exclude it, warn, and keep going

        Returns
        -------
        A list of Classifications instances
        """

        classifications = []

        for classification in self._classifications_from_res_list:
            if self._kwargs["skip_missing"] and not classification.success:
                warnings.warn(
                    "Classification {} not successful. Skipping.".format(classification.id)
                )
                continue

            classifications.append(classification)

        job_names = set([obj.job.name for obj in classifications])

        # warn if some of the classifications in this collection are not alike
        if len(job_names) > 1:
            warnings.warn(
                "SampleCollection contains multiple analysis types: {}".format(", ".join(job_names))
            )

        return classifications

    @cached_property
    def metadata(self) -> pd.DataFrame:
        """Transform a list of Samples or Classifications into a `pd.DataFrame` of metadata."""
        import pandas as pd

        from onecodex.models import Classifications

        DEFAULT_FIELDS = None
        metadata = []

        for obj in self._res_list:
            try:
                classification_id = (
                    obj.id if isinstance(obj, Classifications) else obj.primary_classification.id
                )
            except AttributeError:
                classification_id = None
            sample = obj.sample if isinstance(obj, Classifications) else obj

            m = sample.metadata

            if DEFAULT_FIELDS is None:
                DEFAULT_FIELDS = list(m.__class__.model_fields.keys())
                DEFAULT_FIELDS.remove("field_uri")
                DEFAULT_FIELDS.remove("sample")
                DEFAULT_FIELDS.remove("custom")

            metadatum = {f: getattr(m, f) for f in DEFAULT_FIELDS}
            metadatum["classification_id"] = classification_id
            metadatum["sample_id"] = sample.id
            metadatum["metadata_id"] = m.id
            metadatum["created_at"] = sample.created_at
            metadatum["filename"] = sample.filename
            metadatum["project"] = getattr(sample.project, "name", "")

            metadatum.update(m.custom)
            metadata.append(metadatum)

        if metadata:
            df = pd.DataFrame(metadata)
            index = "classification_id" if df["classification_id"].is_unique else "sample_id"
            metadata = df.set_index(index)
        else:
            metadata = pd.DataFrame(
                columns=["classification_id", "sample_id", "metadata_id", "created_at"]
            )

        return metadata

    def _collate_results(
        self,
        metric: str | Metric = Metric.Auto,
        include_host: bool = False,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Transform a list of Classifications into `pd.DataFrames` of taxonomy and results data.

        Parameters
        ----------
        metric : Metric | str (default=Metric.Auto)
            Which metric to use for the abundance/count of a particular taxon in a sample.
        include_host : bool, optional
            If True, keep (rather than drop) count/abundance data for host taxa (default=False)

        Returns
        -------
        pd.DataFrame of metric values, and a pd.DataFrame of tax_info
        """

        import numpy as np
        import pandas as pd

        metric = Metric.from_value(metric)

        if metric == Metric.Auto:
            metric = self.automatic_metric

        # getting classification IDs is 15% of execution time
        classification_ids = [c.id for c in self._classifications]

        # Compile info about all taxa observed in the classification results
        tax_info = {"tax_id": [], "name": [], "rank": [], "parent_tax_id": []}
        tax_ids = set()

        for classification in self._classifications:
            # pulling results from API is the slowest part of the function, 75% of the execution time
            results = classification.results()
            host_tax_ids = results.get("host_tax_ids", [])

            # d contains info about a taxon in result, including name, id, counts, rank, etc.
            for data in results["table"]:
                d_tax_id = data["tax_id"]

                if not include_host and d_tax_id in host_tax_ids:
                    continue

                if d_tax_id not in tax_ids:
                    for k in ("tax_id", "name", "rank", "parent_tax_id"):
                        tax_info[k].append(data[k])
                    tax_ids.add(d_tax_id)

        tax_info = pd.DataFrame(tax_info, dtype=object, copy=False)
        tax_info.set_index("tax_id", inplace=True)

        # Now that we have a complete list of taxa in these classification results, take a second
        # pass through the results, filling in the abundances/counts into a single numpy array.
        #
        # Note: this is *much* faster than storing the data in individual taxon vectors and letting
        # pandas combine/coerce the data when constructing the dataframe below.
        tax_id_to_idx = {}
        for idx, tax_id in enumerate(tax_info.index):
            tax_id_to_idx[tax_id] = idx

        data = np.zeros((len(self._classifications), len(tax_info.index)), dtype=metric.dtype)

        if AbundanceMetric.has_value(metric):
            data.fill(np.nan)

        for c_idx, c in enumerate(self._classifications):
            # results are cached from the call earlier in this method
            results = c.results()
            host_tax_ids = results.get("host_tax_ids", [])

            for d in results["table"]:
                d_tax_id = d["tax_id"]

                if not include_host and d_tax_id in host_tax_ids:
                    continue

                if metric in (Metric.AbundanceWChildren, Metric.Abundance):
                    data[c_idx, tax_id_to_idx[d_tax_id]] = d.get(metric.results_key)
                else:
                    data[c_idx, tax_id_to_idx[d_tax_id]] = d[metric.results_key] or 0

        df = pd.DataFrame(
            data,
            index=pd.Index(classification_ids, name="classification_id"),
            columns=tax_info.index,
            copy=False,
        )

        return df, tax_info

    @cached_property
    def is_metagenomic(self) -> bool:
        """Return True if this collection appears to be from WGS metagenomics analyses."""

        # warn if some of the classifications in this collection are not alike
        job_names = set([obj.job.name for obj in self._classifications])

        if len(job_names) > 1:
            warnings.warn(
                "SampleCollection contains multiple analysis types: {}".format(", ".join(job_names))
            )

        if len(job_names) == 1 and "One Codex Database" in list(job_names)[0]:
            return True
        else:
            return False

    def automatic_rank(self, metric: str | Metric) -> Rank:
        """Automatically determine the best rank to use for stats/viz.

        If metric is provided and AbundanceWChildren or Abundance -> Species.

        Otherwise, check the job type. If it's metagenomic (as in shotgun) return Species.
        Otherwise, which is probably the case for Targeted Loci analyses, return Genus
        """

        if isinstance(metric, str):
            metric = Metric.from_value(metric)

        if AbundanceMetric.has_value(metric) or self.is_metagenomic:
            return Rank.Species
        else:
            return Rank.Genus

    @cached_property
    def automatic_metric(self) -> Metric:
        """
        Return the best classification metric given the set of classifications.

        Automatically determines the best metric to use. If more than half of the samples have
        abundance estimates, use AbundanceWChildren. Otherwise, use ReadcountWChildren.
        """

        metric_counts = {"abundance": 0, "readcount": 0, "readcount_w_children": 0}

        for classification in self._classifications:
            for row in classification.results()["table"]:
                if row["rank"] != "species":
                    continue
                for key in metric_counts:
                    if key in row and row[key] is not None:
                        metric_counts[key] += 1

        # if abundances are available in more than half of the samples, use that
        if metric_counts["abundance"] >= (metric_counts["readcount_w_children"] // 2):
            return Metric.AbundanceWChildren
        else:
            return Metric.NormalizedReadcountWChildren

    @cached_property
    def _classification_ids_without_abundances(self) -> list[str]:
        """Return classification IDs that have no abundance or abundance_w_children values."""

        classification_ids_without_abundances = []
        for classification in self._classifications:
            has_abundances = False
            for row in classification.results()["table"]:
                if row.get("abundance_w_children") is not None or row.get("abundance") is not None:
                    has_abundances = True
                    break
            if not has_abundances:
                classification_ids_without_abundances.append(classification.id)
        return classification_ids_without_abundances

    @cached_property
    def taxonomy(self):
        # metric is arbitrary here as long as it has a value
        _, tax_info = self._collate_results(metric=Metric.Readcount)

        return tax_info

    @cached_property
    def _functional_profiles(self) -> list[FunctionalProfiles]:
        """Transform a list of Samples or Classifications into a list of FunctionalProfiles objects.

        Each sample will be mapped to the newest job version which may result in mixing different
        result versions on the list.

        Returns
        -------
        A list of FunctionalAnnotations
        """
        from onecodex.models import FunctionalProfiles, Samples

        sample_ids = [
            obj.id if isinstance(obj, Samples) else obj.sample.id for obj in self._res_list
        ]
        # Get all Functional Profiles for the current sample collection
        batch_size = 50
        profiles = []
        for i in range(0, len(sample_ids), batch_size):
            profiles += FunctionalProfiles.where(sample=sample_ids[i : i + batch_size])
        profiles = [fp for fp in profiles if fp.success]

        if not profiles:
            if self._kwargs["skip_missing"]:
                warnings.warn("No functional profiles found for sample collection")
                return []
            else:
                raise OneCodexException("No functional profiles found for sample collection")

        # Determine the newest Functional Analysis version for each sample
        sample_id_to_profile = {}
        for profile in profiles:
            if profile.sample.id in sample_id_to_profile:
                other = sample_id_to_profile[profile.sample.id]
                if other.job.created_at < profile.job.created_at or (
                    other.job.id == profile.job.id and other.created_at < profile.created_at
                ):
                    # Replace if either the job version is older or the run is older
                    sample_id_to_profile[profile.sample.id] = profile
            else:
                sample_id_to_profile[profile.sample.id] = profile

        # Issue missing results or mixed versions warnings
        newest_profiles = list(sample_id_to_profile.values())
        functional_sample_ids = {fp.sample.id for fp in newest_profiles}

        job_ids = {fp.job.id for fp in newest_profiles}
        if len(job_ids) > 1:
            warnings.warn("Be advised: mixing functional profile versions")

        for sample_id in sample_ids:
            if sample_id not in functional_sample_ids:
                msg = f"Functional profile not found for sample {sample_id}."
                if self._kwargs["skip_missing"]:
                    warnings.warn(msg + " Skipping.")
                else:
                    raise OneCodexException(msg)

        return newest_profiles

    @lru_cache
    def _functional_results(
        self,
        annotation: FunctionalAnnotations,
        metric: FunctionalAnnotationsMetric,
        taxa_stratified: bool,
        fill_missing: bool,
        filler: Any,
    ):
        """
        Return a dataframe of all functional profile data and feature id to name mapping.

        Parameters
        ----------
        annotation : {onecodex.lib.enum.FunctionalAnnotations, str}
            Annotation data to return
        taxa_stratified : bool, optional
            Return taxonomically stratified data
        metric : {onecodex.lib.enum.FunctionalAnnotationsMetric, str}
            Metric values to return,
            `{'coverage', 'abundance'}` for `annotation==onecodex.lib.enum.FunctionalAnnotation.Pathways`,
            `{'rpk', 'cpm'}` for other `annotation`
        fill_missing : bool, optional
            Fill `np.nan` values
        filler : float, optional
            Value with which to fill `np.nan` values
        """
        import numpy as np
        import pandas as pd

        # validate args
        annotation = FunctionalAnnotations(annotation)
        metric = FunctionalAnnotationsMetric(metric)

        if annotation == FunctionalAnnotations.Pathways:
            if metric.value not in ["coverage", "abundance", "complete_abundance"]:
                raise ValueError(
                    "If using annotation='pathways', 'metric' must be one of "
                    "['coverage', 'abundance', 'complete_abundance']"
                )
        elif metric.value not in ["cpm", "rpk"]:
            raise ValueError(
                f"If using annotation={annotation.value}, 'metric' must be one of ['cpm', 'rpk']"
            )

        # iterate over functional profiles, subset data, and store in data dict
        data = {}
        all_features = set()
        feature_id_to_name = {}

        sample_id_to_profile_id = {}

        for profile in self._functional_profiles:
            sample_id = profile.sample.id

            if sample_id in sample_id_to_profile_id:
                raise ValueError(f"More than one functional profile for sample {sample_id}")

            sample_id_to_profile_id[sample_id] = profile.id

            # get table using One Codex API
            table = profile.filtered_table(
                annotation=annotation, metric=metric, taxa_stratified=taxa_stratified
            )

            # store tables for later retrieval
            data[sample_id] = dict(zip(table.id, table.value))
            feature_id_to_name.update(dict(zip(table.id, table.name)))
            all_features.update(set(table["id"]))

        features_to_ix = {}
        feature_list = []
        for ix, feature in enumerate(all_features):
            features_to_ix[feature] = ix
            feature_list.append(feature)

        # initialize an array and fill it
        array = np.full(shape=(len(data), len(features_to_ix)), dtype=float, fill_value=np.nan)
        sample_ids = []
        for sample_index, sample_id in enumerate(data):
            for feature_id in data[sample_id]:
                array[sample_index, features_to_ix[feature_id]] = data[sample_id][feature_id]
            sample_ids.append(sample_id)

        functional_profile_ids = [sample_id_to_profile_id[sample_id] for sample_id in sample_ids]

        df = pd.DataFrame(
            array,
            index=pd.Index(functional_profile_ids, name="functional_profile_id"),
            columns=feature_list,
            copy=False,
        )

        if fill_missing:
            df.fillna(filler, inplace=True)

        return df, feature_id_to_name

    def to_otu(
        self,
        biom_id: str | None = None,
        include_ranks: list[str] = CANONICAL_RANKS,
        metric: Metric = Metric.Auto,
    ):
        """
        Generate a BIOM-formatted data structure.

        Parameters
        ----------
        biom_id : string, optional
            Optionally specify an `id` field for the generated v1 BIOM file.

        include_ranks : list
            A list of ranks to include in the taxonomy/OTU table. Uses
            onecodex.models.collection.CANONICAL_RANKS by default.

        Returns
        -------
        otu_table : OrderedDict
            A BIOM OTU table, returned as a Python `OrderedDict` (can be dumped to JSON)
        """
        otu_format = "Biological Observation Matrix 1.0.0"

        # Note: This is exact format URL is required by https://github.com/biocore/biom-format
        otu_url = "http://biom-format.org"

        otu = OrderedDict(
            {
                "id": biom_id,
                "format": otu_format,
                "format_url": otu_url,
                "type": "OTU table",
                "generated_by": "One Codex API V1",
                "date": datetime.now().isoformat(),
                "rows": [],
                "columns": [],
                "matrix_type": "sparse",
                "matrix_element_type": "int",
            }
        )

        rows = defaultdict(dict)

        self._collate_results(include_host=True, metric=metric)  # make sure 9606 is in the taxonomy

        root = self.tree_build()
        for classification in self._classifications:
            col_id = len(otu["columns"])  # 0 index

            # Re-encoding the JSON is a bit of a hack, but
            # we need a ._to_dict() method that properly
            # resolves references and don't have one at the moment
            columns_entry = {
                "id": str(classification.id),
                "sample_id": str(classification.sample.id),
                "sample_filename": classification.sample.filename,
                "metadata": json.loads(
                    classification.sample.metadata.model_dump_json(exclude={"field_uri", "sample"})
                ),
            }

            otu["columns"].append(columns_entry)
            sample_df = classification.table()

            for row in sample_df.iterrows():
                tax_id = row[1]["tax_id"]
                tax_node = root.find(tax_id)

                # only keep canonical rows (e.g., don't include things like
                # "root" in the OTU table)
                if tax_node.rank not in include_ranks:
                    continue

                rows[tax_id][col_id] = int(row[1][Metric.Readcount])

        num_rows = len(rows)
        num_cols = len(otu["columns"])

        otu["shape"] = [num_rows, num_cols]
        otu["data"] = []

        for tax_id in sorted(rows):
            row_id = len(otu["rows"])
            tax_node = root.find(tax_id)
            lineage = tax_node.ancestors() + [tax_node]

            rank_to_name = {n.rank: n.tax_name for n in lineage}

            canonical_names = []
            for rank in CANONICAL_RANKS:
                canonical_names.append(rank_to_name.get(rank, ""))

            otu["rows"].append({"id": tax_id, "metadata": {"taxonomy": canonical_names}})

            for sample_with_hit in rows[tax_id]:
                counts = rows[tax_id][sample_with_hit]
                otu["data"].append([row_id, sample_with_hit, counts])

        return otu

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
        rank: Rank | str = Rank.Auto,
        top_n: Optional[int] = None,
        threshold: Optional[float] = None,
        remove_zeros: bool = True,
        include_host: bool = False,
        table_format: Literal["wide", "long"] = "wide",
        include_taxa_missing_rank: bool = False,
        fill_missing: bool = True,
        filler: Any = 0,
        metric: Metric = Metric.Auto,
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
        metric: Metric
            the taxonomic abundance metric to use. See onecodex.lib.enums.Metric for definitions
        threshold : `float`, optional
            Return only taxa more abundant than this threshold in one or more samples.
        remove_zeros : `bool`, optional
            Do not return taxa that have zero abundance in every sample.
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
        return self._to_classification_df(
            rank=rank,
            metric=metric,
            top_n=top_n,
            threshold=threshold,
            remove_zeros=remove_zeros,
            table_format=table_format,
            include_taxa_missing_rank=include_taxa_missing_rank,
            fill_missing=fill_missing,
            include_host=include_host,
            filler=filler,
        ).copy()

    @lru_cache
    def _to_classification_df(
        self,
        rank: Rank | str = Rank.Auto,
        metric: Metric | str = Metric.Auto,
        include_host: bool = False,
        top_n: Optional[int] = None,
        threshold: Optional[float] = None,
        remove_zeros: bool = True,
        table_format: Literal["wide", "long"] = "wide",
        include_taxa_missing_rank: bool = False,
        fill_missing: bool = True,
        filler: Any = 0,
    ) -> "ClassificationsDataFrame":
        """Generate the ClassificationsDataFrame from data from the One Codex API."""

        from onecodex.dataframes import ClassificationsDataFrame

        metric, rank = self._parse_classification_config_args(metric=metric, rank=rank)

        if include_taxa_missing_rank:
            if rank is None:
                raise OneCodexException(
                    "`rank` must be specified when passing `include_taxa_missing_rank=True`."
                )

        # after this point, _to_classification_df should no longer care about metric or include_host
        df, tax_info = self._collate_results(
            metric=metric,
            include_host=include_host,
        )

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
                    tax_id_level = Rank(tax_info["rank"][tax_id]).level
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

        if metric in (
            Metric.PropReadcount,
            Metric.PropReadcountWChildren,
            Metric.PropClassified,
            Metric.PropClassifiedWChildren,
            Metric.NormalizedReadcount,
            Metric.NormalizedReadcountWChildren,
        ):
            if metric in (Metric.PropClassified, Metric.PropClassifiedWChildren):
                denoms = [
                    c._classification_stats["n_mapped_microbial_reads"]
                    for c in self._classifications
                ]
            elif metric in (Metric.PropReadcount, Metric.PropReadcountWChildren):
                denoms = [c._classification_stats["n_reads_total"] for c in self._classifications]
            elif metric in (Metric.NormalizedReadcount, Metric.NormalizedReadcountWChildren):
                denoms = df.sum(axis=1)
            else:
                raise Exception("unreachable")

            df = df.div(denoms, axis=0).fillna(0.0)

        if fill_missing:
            df = df.fillna(filler)

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
            "ocx_metric": metric,
            # this does not include parents...
            "ocx_taxonomy": tax_info.copy(),
            "ocx_threshold": threshold,
            "ocx_classification_ids_without_abundances": self._classification_ids_without_abundances,
        }

        # generate long-format table
        if table_format == "long":
            metric_display_name = self._display_name_for_metric(metric=metric)

            long_df = {"classification_id": [], "tax_id": [], metric_display_name: []}

            for t_id in df:
                for c_id, count in df[t_id].items():
                    long_df["classification_id"].append(c_id)
                    long_df["tax_id"].append(t_id)
                    long_df[metric_display_name].append(count)

            results_df = ClassificationsDataFrame(long_df, **ocx_data)
        elif table_format == "wide":
            results_df = ClassificationsDataFrame(df, **ocx_data)
        else:
            raise OneCodexException("table_format must be one of: long, wide")

        return results_df

    def _metadata_fetch(
        self,
        metadata_fields,
        results_df: pd.DataFrame,
        label=None,
        match_taxonomy=True,
        coerce_missing_composite_fields=True,
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

        magic_metadata = pd.DataFrame(
            {"classification_id": [c.id for c in self._classifications]}
        ).set_index("classification_id")

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
                elif match_taxonomy and str_f in results_df.keys():
                    # is a tax_id
                    rank = self.taxonomy["rank"][str_f]
                    if Rank.has_value(rank):
                        tax_name = self.taxonomy["name"][str_f]
                        renamed_field = "{} ({})".format(tax_name, str_f)
                        magic_metadata[renamed_field] = results_df[str_f]
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
                        tax_id, name = hits[0]
                        # report within-rank abundance (we have to fetch this from a new
                        # results_df because it may be at a different rank)
                        df = self.to_df(
                            metric=results_df.ocx_metric,
                            rank=self.taxonomy.loc[tax_id]["rank"],
                        )

                        renamed_field = f"{name} ({tax_id})"
                        magic_metadata[renamed_field] = df[tax_id]
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

    def _display_name_for_metric(self, metric: Metric) -> str:
        if metric == Metric.Auto:
            metric = self.automatic_metric
        return metric.display_name
