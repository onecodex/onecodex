from __future__ import annotations

import warnings
from typing import Callable, Any, TYPE_CHECKING

from onecodex.lib.enums import (
    Metric,
    Rank,
    FunctionalAnnotations,
    FunctionalAnnotationsMetric,
    Link,
)
from onecodex.models import SampleCollection as BaseSampleCollection
from onecodex.exceptions import (
    OneCodexException,
    PlottingException,
    PlottingWarning,
    ValidationError,
)
from .enums import PlotRepr, PlotType
from .metadata import metadata_record_to_label, deduplicate_labels, sort_metadata_records
from .utils import get_plot_title
from .models import PlotParams, PlotResult
from .export import export_chart_data

if TYPE_CHECKING:
    import pandas as pd


METADATA_FIELD_PLOT_PARAMS = [
    "facet_by",
    "group_by",
    "secondary_group_by",
    "filter_by",
    "label_by",
    "sort_by",
]


###
# SampleCollection shims to support data fetched via the internal v2 API instead of the v1 API. The
# v2 API has an endpoint that provides all the necessary/minimal sample data, and is more efficient
# than fetching the data via the v1 API.
###


class Samples:
    """Mock the Samples model."""

    def __init__(self, sample_datum: dict):
        self._sample_datum = sample_datum

    @property
    def id(self) -> str:
        return self._sample_datum["uuid"]

    @property
    def _metadata(self) -> dict:
        return self._sample_datum["metadata"]

    @property
    def _classification(self) -> dict | None:
        return self._sample_datum["primary_classification"]

    @property
    def _functional_profile(self) -> dict | None:
        return self._sample_datum["functional_profile"]


class Classifications(dict):
    """Mock the Classifications model."""

    id: str | None = None

    def results(self) -> "Classifications":
        return self


class FunctionalProfiles:
    """Mock the FunctionalProfiles model."""

    def __init__(self, functional_run_uuid: str, sample_uuid: str, results: dict):
        self._uuid = functional_run_uuid
        self._results = results
        self._sample_uuid = sample_uuid

    @property
    def id(self) -> str:
        return self._uuid

    @property
    def sample(self) -> Samples:
        return Samples({"uuid": self._sample_uuid})

    def _filtered_results(
        self,
        annotation: FunctionalAnnotations,
        metric: FunctionalAnnotationsMetric,
        taxa_stratified: bool,
    ) -> dict:
        if taxa_stratified:
            raise OneCodexException("Taxa stratified results are not currently supported")

        return {
            "table": self._results.get(f"{annotation}-{metric}", []),
            "n_reads": self._results["n_reads"],
            "n_mapped": self._results["n_mapped"],
        }

    def filtered_table(self, *args, **kwargs) -> pd.DataFrame:
        from onecodex.models import FunctionalProfiles

        return FunctionalProfiles.filtered_table(self, *args, **kwargs)


class SampleCollection(BaseSampleCollection):
    def __init__(self, samples: list[Samples], **kwargs):
        """Overridden for shims."""
        # For some reason, it would try to collate_results in loop because that value is not set.
        # Strangely, that value should be set in collate_results so it shouldn't recurse.
        # Not sure what's happening but setting it initially fixes it.
        self._cached = {"metric": kwargs.get("metric", Metric.Auto)}
        self._kwargs = {
            "skip_missing": True,
            "metric": Metric.Auto,
            "include_host": False,
            "job": None,
        }
        self._kwargs.update(kwargs)

        self.samples = samples
        self._res_list = self.samples
        self._collate_metadata()

        self._classification_fetch()
        self.rank = self._default_rank
        self._functional_profiles  # precompute the cache

    def _collate_metadata(self):
        """Overridden for shims."""
        import pandas as pd

        metadata = [sample._metadata for sample in self.samples]

        if metadata:
            df = pd.DataFrame(metadata)
            index = "classification_id" if df["classification_id"].is_unique else "sample_id"
            metadata = df.set_index(index)
        else:
            metadata = pd.DataFrame(
                columns=["classification_id", "sample_id", "metadata_id", "created_at"]
            )

        self._cached["metadata"] = metadata

    def _classification_fetch(self):
        """Overridden for shims."""
        classifications = []
        job_names = set()
        for sample in self.samples:
            summary = sample._classification
            if not summary:
                continue

            results = Classifications()
            results.update(summary["api_results"])
            results["id"] = summary["uuid"]
            results.id = summary["uuid"]

            job_names.add(summary["job_name"])
            classifications.append(results)

        self._cached["is_metagenomic"] = False
        if all(("One Codex Database" in name for name in job_names)):
            self._cached["is_metagenomic"] = True

        self._default_rank = self._get_auto_rank(Rank.Auto)
        self._cached["classifications"] = classifications

    @property
    def _functional_profiles(self):
        """Overridden for shims."""
        if "functional_profiles" not in self._cached:
            functional_results = []
            for sample in self.samples:
                profile = sample._functional_profile
                if not profile:
                    continue
                functional_results.append(
                    FunctionalProfiles(profile["uuid"], profile["sample_uuid"], profile["results"])
                )

            self._cached["functional_profiles"] = functional_results

        return self._cached["functional_profiles"]

    def plot(self, params: PlotParams) -> PlotResult:
        import altair as alt

        result = None
        with warnings.catch_warnings(record=True) as captured_warnings:
            warnings.simplefilter("always", PlottingWarning)

            try:
                result = self._plot(params)
            except (ValidationError, PlottingException) as e:
                # Expected user error
                return PlotResult(params=params, error=str(e))
            except alt.MaxRowsError:
                return PlotResult(
                    params=params,
                    error="The selected dataset is too large to plot. Please try a different plot type or select a fewer number of samples.",
                )

        for warning in captured_warnings:
            if warning.category is PlottingWarning:
                result.warnings.append(str(warning.message))
            else:
                warnings.warn(warning.message, warning.category)

        return result

    def _plot(self, params: PlotParams) -> PlotResult:
        self._validate_plot_params(params)

        if params.filter_by and params.filter_value:
            # Create a *new* filtered SampleCollection and reassign to `self`, rather than filtering
            # the current `self` in-place. We don't want to cache the filtered SampleCollection
            self = self._filter_by_metadata(params.filter_by, params.filter_value)

        if params.metric == Metric.Auto:
            self._collate_results()  # This resolves `self._metric`
            params = params.model_copy(
                update={"metric": Metric(self._metric)}
            )  # don't mutate the input

        label_func = self._x_axis_label_func(params.plot_type, params.label_by)
        if params.plot_type == PlotType.Functional:
            x_axis_label_links = self._x_axis_label_functional_links(label_func, params.group_by)
        else:
            x_axis_label_links = self._x_axis_label_classification_links(
                label_func, params.group_by
            )

        sort_x_func = self._x_axis_sort_func(params.sort_by, label_func)
        title = get_plot_title(params)
        default_x_axis_title = "Samples"
        # "container" for responsive plots when window is resized
        default_size_kwargs = {"width": "container", "height": "container"}

        if params.plot_type == PlotType.Taxa:
            if params.facet_by:
                # "container" doesn't currently work with facet plots in vega-lite/altair, so fall
                # back to default size (DEV-4753)
                default_size_kwargs = {}

            if params.plot_repr == PlotRepr.Bargraph:
                chart = self.plot_bargraph(
                    return_chart=True,
                    top_n=params.top_n,
                    rank=params.rank,
                    haxis=params.facet_by,
                    title=title,
                    xlabel=None if params.facet_by or params.group_by else default_x_axis_title,
                    label=None if params.group_by else label_func,
                    sort_x=None if params.group_by else sort_x_func,
                    group_by=params.group_by,
                    link=Link.Ncbi,
                    match_taxonomy=False,
                    **default_size_kwargs,
                )
            else:
                chart = self.plot_heatmap(
                    return_chart=True,
                    top_n=params.top_n,
                    rank=params.rank,
                    haxis=params.facet_by,
                    title=title,
                    xlabel=None if params.facet_by else default_x_axis_title,
                    label=label_func,
                    sort_x=sort_x_func,
                    link=Link.Ncbi,
                    match_taxonomy=False,
                    **default_size_kwargs,
                )
        elif params.plot_type == PlotType.Alpha:
            if params.facet_by:
                # "container" doesn't currently work with facet plots in vega-lite/altair, so fall
                # back to default size (DEV-4753)
                default_size_kwargs = {}

            if params.facet_by:
                xlabel = None
            elif params.group_by:
                xlabel = params.group_by
            else:
                xlabel = default_x_axis_title

            chart = self.plot_metadata(
                return_chart=True,
                rank=params.rank,
                vaxis=params.alpha_metric,
                haxis=params.group_by or "Label",
                secondary_haxis=params.secondary_group_by,
                facet_by=params.facet_by,
                title=title,
                xlabel=xlabel,
                label=label_func,
                sort_x=sort_x_func,
                coerce_haxis_dates=False,  # dates formatted by Custom Plots look nicer
                match_taxonomy=False,
                **default_size_kwargs,
            )
        elif params.plot_type == PlotType.Beta:
            if params.plot_repr == PlotRepr.Pcoa:
                chart = self.plot_mds(
                    return_chart=True,
                    rank=params.rank,
                    metric=params.beta_metric,
                    color=params.facet_by,
                    title=title,
                    label=label_func,
                    match_taxonomy=False,
                    **default_size_kwargs,
                )
            elif params.plot_repr == PlotRepr.Pca:
                chart = self.plot_pca(
                    return_chart=True,
                    rank=params.rank,
                    color=params.facet_by,
                    title=title,
                    label=label_func,
                    match_taxonomy=False,
                    **default_size_kwargs,
                )
            elif params.plot_repr == PlotRepr.Distance:
                # "container" doesn't currently work with compound plots in vega-lite/altair, so
                # fall back to default size (DEV-4753)
                default_size_kwargs = {}
                chart = self.plot_distance(
                    return_chart=True,
                    rank=params.rank,
                    metric=params.beta_metric,
                    title=title,
                    xlabel=default_x_axis_title,
                    label=label_func,
                    match_taxonomy=False,
                    **default_size_kwargs,
                )

        elif params.plot_type == PlotType.Functional:
            if params.facet_by:
                # "container" doesn't currently work with facet plots in vega-lite/altair, so fall
                # back to default size (DEV-4753)
                default_size_kwargs = {}
            else:
                if params.functional_top_n > 30:
                    default_size_kwargs = {"width": "container"}

            functional_metric = params.functional_metric
            if params.functional_annotation == FunctionalAnnotations.Pathways:
                functional_metric = params.functional_pathways_metric
            chart = self.plot_functional_heatmap(
                return_chart=True,
                title=title,
                annotation=params.functional_annotation,
                metric=functional_metric,
                top_n=params.functional_top_n,
                sort_x=sort_x_func,
                label=label_func,
                function_label=params.functional_label,
                haxis=params.facet_by,
                xlabel=None if params.facet_by else default_x_axis_title,
                **default_size_kwargs,
            )
        else:
            raise OneCodexException(f"Unknown plot type: {params.plot_type}")

        # Open links in new tab: https://stackoverflow.com/a/72241020/3776794
        chart["usermeta"] = {"embedOptions": {"loader": {"target": "_blank", "rel": "noreferrer"}}}

        exported_chart_data = export_chart_data(params, chart)

        # This is a backwards compatibility fix.
        # Default OCX plot styles include background and no grid. Custom Plots historically
        # had no background and enabled grid which was due to an unexpected error in loading
        # the `altair` module. This function removes some of the default styling to keep the
        # charts consistent.
        chart = chart.to_dict()
        if isinstance(chart.get("config"), dict):
            chart["config"].pop("background", None)
            if isinstance(chart["config"].get("axis"), dict):
                chart["config"]["axis"].pop("grid", None)

        return PlotResult(
            params=params,
            chart=chart,
            x_axis_label_links=x_axis_label_links,
            exported_chart_data=exported_chart_data,
        )

    def _validate_plot_params(self, params: PlotParams):
        if params.plot_type == PlotType.Functional:
            if not self._functional_profiles:
                raise ValidationError(
                    "Functional Analysis has not been run for any of the selected samples."
                )
        else:
            if not self._classifications:
                raise ValidationError(
                    "Classification results are not available for any of the selected samples."
                )

        for attr in METADATA_FIELD_PLOT_PARAMS:
            fields = getattr(params, attr)
            if not isinstance(fields, list):
                fields = [fields]
            for field in fields:
                if field is not None and field not in self.metadata.columns:
                    attr_display_name = attr.replace("_", " ").title()
                    raise ValidationError(
                        f"The metadata field {field!r} does not exist. Please select a valid "
                        f"metadata field in the {attr_display_name} dropdown."
                    )

    def _filter_by_metadata(self, field: str, values_to_keep: list[str]) -> "SampleCollection":
        import pandas as pd

        def _filter_func(sample: Samples) -> bool:
            metadata = self.metadata
            if field not in metadata.columns:
                return False
            if metadata.index.name == "sample_id":
                metadatum = metadata.loc[sample.id, field]
            else:
                metadatum = metadata.loc[metadata["sample_id"] == sample.id, field].iloc[0]
            return not pd.isna(metadatum) and metadatum in values_to_keep

        return self.filter(_filter_func)

    def _x_axis_label_func(self, plot_type: PlotType, label_by: list[str]) -> Callable[[dict], str]:
        id_field = None
        ids = set()
        if plot_type == PlotType.Functional:
            id_field = "sample_id"
            for profile in self._functional_profiles:
                ids.add(profile.sample.id)
        else:
            id_field = "classification_id"
            for classification in self._classifications:
                ids.add(classification.id)

        labels_by_metadata_id = {}
        for idx, record in self.metadata.to_dict("index").items():
            id_ = record[id_field] if self.metadata.index.name != id_field else idx
            if id_ in ids:
                # Only deduplicate labels for samples that will be included in the plot: if it's a
                # functional plot, only deduplicate labels of samples that have functional profiles.
                # If it's not a functional plot, only deduplicate labels of samples that have
                # classifications.
                labels_by_metadata_id[record["metadata_id"]] = metadata_record_to_label(
                    record, label_by
                ).strip()

        unique_labels_by_metadata_id = deduplicate_labels(labels_by_metadata_id)

        def _label_func(record: dict) -> str:
            return unique_labels_by_metadata_id.get(record["metadata_id"], "N/A")

        return _label_func

    def _x_axis_label_functional_links(
        self, label_func: Callable[[dict], str], group_by: str | None
    ) -> dict[str, str]:
        if group_by:
            return {}

        x_axis_label_links = {}
        sample_uuid_to_functional_uuid = {
            profile.sample.id: profile.id for profile in self._functional_profiles
        }

        # Map each x-axis label to a link containing its functional analysis results.
        for idx, record in self.metadata.to_dict("index").items():
            sample_uuid = record["sample_id"] if self.metadata.index.name != "sample_id" else idx
            if sample_uuid not in sample_uuid_to_functional_uuid:
                continue
            label = label_func(record)
            assert label not in x_axis_label_links
            x_axis_label_links[label] = f"/functional/{sample_uuid_to_functional_uuid[sample_uuid]}"

        return x_axis_label_links

    def _x_axis_label_classification_links(
        self, label_func: Callable[[dict], str], group_by: str | None
    ) -> dict[str, str]:
        if group_by:
            return {}

        # Map each x-axis label to a link containing its classification results.
        x_axis_label_links = {}
        for classification_id, record in self.metadata.to_dict("index").items():
            label = label_func(record)
            assert label not in x_axis_label_links
            x_axis_label_links[label] = f"/classification/{classification_id}"

        return x_axis_label_links

    def _x_axis_sort_func(
        self, sort_by: str | None, label_func: Callable[[dict], str]
    ) -> Callable[[Any], list[str]] | None:
        if sort_by is None:
            return None

        def _sort_x_func(_: Any) -> list[str]:
            records = self.metadata.to_dict("records")
            sorted_records = sort_metadata_records(records, sort_by)
            return [label_func(x) for x in sorted_records]

        return _sort_x_func
