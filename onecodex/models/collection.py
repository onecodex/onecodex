from collections import defaultdict, OrderedDict
from collections.abc import MutableSequence
from typing import Any
from datetime import datetime
import json
import warnings

from onecodex.exceptions import OneCodexException
from onecodex.lib.enums import (
    AbundanceMetric,
    Metric,
    FunctionalAnnotations,
    FunctionalAnnotationsMetric,
)

try:
    from onecodex.analyses import AnalysisMixin
except ImportError:

    class AnalysisMixin(object):
        pass


from onecodex.models.sample import Samples
from onecodex.models.analysis import Classifications

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


class SampleCollection(AnalysisMixin, MutableSequence):
    """A collection of `Samples` or `Classifications` objects.

    Includes lots of methods for analysis of classifications results.

    Notes
    -----
    Inherits from `ResourceList` to provide a list-like API, and `AnalysisMixin` to provide relevant
    analysis methods.
    """

    def __init__(self, *args, **kwargs):
        """Instantiate a new SampleCollection containing `Samples` or `Classifications` objects.

        Parameters
        ----------
        objects : list
            A list of `onecodex.models.Samples` or `onecodex.models.Classifications` objects
            which will be processed into a `SampleCollection`

        skip_missing : bool, optional
            If an analysis was not successful, exclude it, warn, and keep going

        metric : {'readcount_w_children', 'readcount', 'abundance'}, optional
            Which metric to use for the abundance/count of a particular taxon in a sample.

            - 'readcount_w_children': total reads of this taxon and all its descendants
            - 'readcount': total reads of this taxon
            - 'abundance': genome size-normalized relative abundances, from shotgun sequencing

        include_host : bool, optional
            If True, keep (rather than drop) count/abundance data for host taxa

        Examples
        --------
        Given a list of Samples, create a new SampleCollection using abundances:

            samples = [sample1, sample2, sample3]
            collection = SampleCollection(samples, metric='abundance')

        Notes
        -----
        To provide access to the list-like API of `ResourceList`, must also accept a list of
        unwrapped potion resources and a One Codex model.
        """

        if "field" in kwargs and "metric" in kwargs:
            raise OneCodexException(
                "Cannot provide both `field` and `metric`. `field` has been deprecated in favor of "
                "`metric`."
            )
        if "field" in kwargs:
            warnings.warn(
                "The `field` parameter has been renamed to `metric`. Passing `field` to a "
                "SampleCollection is deprecated and will be removed in a future release.",
                DeprecationWarning,
            )
            kwargs["metric"] = kwargs.pop("field")

        self._kwargs = kwargs
        include_host = kwargs.get("include_host", False)
        metric = kwargs.get("metric", "auto")
        job = kwargs.get("job", None)
        skip_missing = kwargs.get("skip_missing", True)
        objects = [*args[0]]
        _optional_model = args[1] if len(args) > 1 else None  # For backwards compatibility

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
            "metric": metric,
            "include_host": include_host,
            "job": job,
        }
        self._cached = {}
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

    def _classification_fetch(self, skip_missing=None):
        """Transform a list of Samples or Classifications into a list of Classifications objects.

        Parameters
        ----------
        skip_missing : bool
            If an analysis was not successful, exclude it, warn, and keep going

        Returns
        -------
        None, but stores a result in self._cached.
        """
        from onecodex.models import Classifications, Samples

        skip_missing = skip_missing if skip_missing else self._kwargs["skip_missing"]

        new_classifications = []

        for obj in self._res_list:
            if isinstance(obj, Samples):
                classification = obj.primary_classification
            elif isinstance(obj, Classifications):
                classification = obj
            else:
                raise OneCodexException(
                    "Objects in SampleCollection must be one of: Classifications, Samples"
                )

            if not classification:
                msg = f"Classification not found for sample {obj.id}."
                if skip_missing:
                    warnings.warn(msg + " Skipping.")
                    continue
                else:
                    raise OneCodexException(msg)

            if skip_missing and not classification.success:
                warnings.warn(f"Classification {classification.id} not successful. Skipping.")
                continue

            new_classifications.append(classification)

        # warn if some of the classifications in this collection are not alike
        job_names = set([obj.job.name for obj in new_classifications])

        if len(job_names) > 1:
            warnings.warn(
                "SampleCollection contains multiple analysis types: {}".format(", ".join(job_names))
            )

        self._cached["is_metagenomic"] = False
        if len(job_names) == 1 and "One Codex Database" in list(job_names)[0]:
            self._cached["is_metagenomic"] = True

        self._cached["classifications"] = new_classifications

    @property
    def _classifications(self):
        if "classifications" not in self._cached:
            self._classification_fetch()

        return self._cached["classifications"]

    def _collate_metadata(self):
        """Transform a list of Samples or Classifications into a `pd.DataFrame` of metadata."""
        from onecodex.models import Classifications
        import pandas as pd

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

            # Standard/structured metadata fields get priority over custom metadata fields in the
            # event of a name clash. If the user is allowed to override certain metadata fields like
            # `sample_id`, `metadata_id`, and `classification_id`, downstream code will break in
            # strange ways.
            m = sample.metadata
            metadatum = {}
            metadatum.update(m.custom)

            if DEFAULT_FIELDS is None:
                DEFAULT_FIELDS = list(m.__class__.model_fields.keys())
                DEFAULT_FIELDS.remove("field_uri")
                DEFAULT_FIELDS.remove("sample")
                DEFAULT_FIELDS.remove("custom")

            for f in DEFAULT_FIELDS:
                metadatum[f] = getattr(m, f)

            metadatum["classification_id"] = classification_id
            metadatum["sample_id"] = sample.id
            metadatum["metadata_id"] = m.id
            metadatum["created_at"] = sample.created_at
            metadatum["filename"] = sample.filename
            metadatum["project"] = getattr(sample.project, "name", "")
            metadatum["sample_name"] = m.name or sample.filename or "N/A"

            metadata.append(metadatum)

        if metadata:
            df = pd.DataFrame(metadata)
            index = "classification_id" if df["classification_id"].is_unique else "sample_id"
            metadata = df.set_index(index)
        else:
            metadata = pd.DataFrame(
                columns=["classification_id", "sample_id", "metadata_id", "created_at"]
            )

        self._cached["metadata"] = metadata

    @property
    def metadata(self):
        if "metadata" not in self._cached:
            self._collate_metadata()

        return self._cached["metadata"]

    def _collate_results(self, metric=None, include_host=None):
        """Transform a list of Classifications into `pd.DataFrames` of taxonomy and results data.

        Parameters
        ----------
        metric : {'readcount_w_children', 'readcount', 'abundance'}
            Which metric to use for the abundance/count of a particular taxon in a sample.

            - 'readcount_w_children': total reads of this taxon and all its descendants
            - 'readcount': total reads of this taxon
            - 'abundance': genome size-normalized relative abundances, from shotgun sequencing

        include_host : bool, optional
            If True, keep (rather than drop) count/abundance data for host taxa

        Returns
        -------
        None, but stores a result in self._cached.
        """
        import numpy as np
        import pandas as pd

        metric = metric if metric else self._kwargs["metric"]
        include_host = include_host if include_host else self._kwargs["include_host"]

        try:
            metric = Metric(metric)
        except ValueError:
            raise OneCodexException("Specified metric ({}) not valid.".format(metric))

        # getting classification IDs is 15% of execution time
        classification_ids = [c.id for c in self._classifications]

        # Compile info about all taxa observed in the classification results & get classification_ids_without_abundances
        tax_info = {"tax_id": [], "name": [], "rank": [], "parent_tax_id": []}
        tax_ids = set()
        classification_ids_without_abundances = []

        for classification in self._classifications:
            # pulling results from API is the slowest part of the function, 75% of the execution time
            results = classification.results()
            host_tax_ids = results.get("host_tax_ids", [])
            has_abundance_estimates = False

            # d contains info about a taxon in result, including name, id, counts, rank, etc.
            for data in results["table"]:
                d_tax_id = data["tax_id"]

                if not include_host and d_tax_id in host_tax_ids:
                    continue

                if (
                    data.get(Metric.AbundanceWChildren) is not None
                    or data.get(Metric.Abundance) is not None
                ):
                    has_abundance_estimates = True

                if d_tax_id not in tax_ids:
                    for k in ("tax_id", "name", "rank", "parent_tax_id"):
                        tax_info[k].append(data[k])
                    tax_ids.add(d_tax_id)

            if not has_abundance_estimates:
                classification_ids_without_abundances.append(classification.id)

        if metric == Metric.Auto:
            metric = Metric.ReadcountWChildren
            # If half or fewer than half of total samples lack abundance estimates, we are ok to
            # use abundance with children as metric. Otherwise, use readcount w children
            if self._is_metagenomic and len(classification_ids_without_abundances) <= (
                len(self) / 2
            ):
                metric = Metric.AbundanceWChildren

        metric_dtype = metric.dtype
        metric = metric.value
        self._cached["metric"] = metric
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

        data = np.zeros((len(self._classifications), len(tax_info.index)), dtype=metric_dtype)

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

                if metric in [Metric.AbundanceWChildren, Metric.Abundance]:
                    data[c_idx, tax_id_to_idx[d_tax_id]] = d[metric]
                else:
                    data[c_idx, tax_id_to_idx[d_tax_id]] = d[metric] or 0

        df = pd.DataFrame(
            data,
            index=pd.Index(classification_ids, name="classification_id"),
            columns=tax_info.index,
            copy=False,
        )

        self._cached["results"] = df
        self._cached["taxonomy"] = tax_info

    @property
    def _metric(self):
        if "metric" not in self._cached:
            self._collate_results()

        return self._cached["metric"]

    @property
    def _is_metagenomic(self):
        if "is_metagenomic" not in self._cached:
            # It looks like ._collate_results() does not actually populate 'is_metagenomic';
            # should call _classification_fetch() instead?
            self._collate_results()

        return self._cached["is_metagenomic"]

    @property
    def _results(self):
        if "results" not in self._cached:
            self._collate_results()

        return self._cached["results"]

    @property
    def taxonomy(self):
        if "taxonomy" not in self._cached:
            self._collate_results()

        return self._cached["taxonomy"]

    def _functional_profiles_fetch(self, skip_missing=None):
        """Transform a list of Samples or Classifications into a list of FunctionalProfiles objects.

        Each sample will be mapped to the newest job version which may result in mixing different
        result versions on the list.

        Parameters
        ----------
        skip_missing : bool
            If an analysis is missing or was not successful, exclude it, warn, and keep going

        Returns
        -------
        None, but stores a result in self._cached.
        """
        from onecodex.models import Samples, FunctionalProfiles

        skip_missing = skip_missing if skip_missing else self._kwargs["skip_missing"]

        sample_ids = [
            obj.id if isinstance(obj, Samples) else obj.sample.id for obj in self._res_list
        ]

        # Get all Functional Profiles for the current sample collection, and determine the newest
        # Functional Analysis version for each sample
        batch_size = 50
        sample_id_to_profile = {}
        for i in range(0, len(sample_ids), batch_size):
            for profile in FunctionalProfiles.where(sample=sample_ids[i : i + batch_size]):
                if skip_missing and not profile.success:
                    warnings.warn(f"Functional profile {profile.id} not successful. Skipping.")
                    continue

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
                if skip_missing:
                    warnings.warn(msg + " Skipping.")
                else:
                    raise OneCodexException(msg)

        self._cached["functional_profiles"] = newest_profiles

    @property
    def _functional_profiles(self):
        if "functional_profiles" not in self._cached:
            self._functional_profiles_fetch()

        return self._cached["functional_profiles"]

    def _collate_functional_results(
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

    def _functional_results(self, **kwargs):
        result_key = "functional_results_" + "_".join([f"{k}={v}" for k, v in kwargs.items()])
        feature_map_key = "functional_features_" + "_".join([f"{k}={v}" for k, v in kwargs.items()])
        if result_key not in self._cached or feature_map_key not in self._cached:
            df, feature_map = self._collate_functional_results(**kwargs)
            self._cached[result_key] = df
            self._cached[feature_map_key] = feature_map
        return self._cached[result_key], self._cached[feature_map_key]

    def to_otu(self, biom_id=None, include_ranks=CANONICAL_RANKS):
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

        self._collate_results(include_host=True)  # make sure 9606 is in the taxonomy

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
