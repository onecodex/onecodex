from __future__ import annotations

import warnings
from functools import cached_property
from typing import Any, Callable

from pydantic import PrivateAttr

from onecodex.exceptions import (
    NoTaxaException,
    OneCodexException,
    OneCodexUserWarning,
    PlottingException,
    StatsException,
    ValidationError,
)
from onecodex.lib.enums import (
    FunctionalAnnotations,
    Link,
    Metric,
)
from onecodex.models import Classifications as BaseClassifications
from onecodex.models import FunctionalProfiles as BaseFunctionalProfiles
from onecodex.models import Jobs, Metadata, Samples
from onecodex.models import SampleCollection as BaseSampleCollection

from .enums import PlotRepr, PlotType, StatsType
from .export import export_chart_data
from .metadata import deduplicate_labels, metadata_record_to_label, sort_metadata_records
from .models import BaseParams, PlotParams, PlotResults, StatsParams, StatsResults
from .utils import get_plot_title

METADATA_FIELD_PLOT_PARAMS = [
    "facet_by",
    "group_by",
    "secondary_group_by",
    "filter_by",
    "label_by",
    "sort_by",
]

METADATA_FIELD_STATS_PARAMS = [
    "group_by",
    "secondary_group_by",
    "filter_by",
    "paired_by",
]


###
# Custom Plots runs in Pyodide. Samples and analyses models are built from the frontend API endpoint
# response (`/api/frontend/custom-plots/sample-data`), and results are downloaded asynchronously
# before plotting.
###


class Classifications(BaseClassifications):
    _loaded_results: dict | None = PrivateAttr(default=None)

    def _results(self) -> dict:
        if self._loaded_results is None:
            raise OneCodexException(f"Results have not been loaded for classification {self.id}")
        return self._loaded_results


class FunctionalProfiles(BaseFunctionalProfiles):
    _loaded_results: dict | None = PrivateAttr(default=None)

    def _condensed_results(self) -> dict | None:
        results = self._loaded_results
        if results is None or results.get("version") != self._FUNCTIONAL_RESULTS_VERSION:
            return None
        return results


def samples_from_sample_data(
    sample_data: list[dict],
) -> tuple[list[Samples], dict[str, FunctionalProfiles]]:
    """Build models from `/api/frontend/custom-plots/sample-data` response.

    In order to skip model validation, `model_construct` is used (not all model fields
    are returned from the endpoint). References are set to model objects as `ApiRef`
    cannot be resolved. All metadata is set as `custom`.

    Returns the samples and their functional profiles keyed by sample ID.
    """
    samples = []
    functional_profiles = {}

    for datum in sample_data:
        metadata = datum["metadata"]
        sample = Samples.model_construct(
            field_uri=Samples._convert_id_to_uri(datum["uuid"]),
            created_at=metadata.get("created_at"),
            metadata=Metadata.model_construct(
                field_uri=Metadata._convert_id_to_uri(metadata["metadata_id"]),
                custom=metadata,
            ),
        )

        summary = datum.get("primary_classification")
        if summary:
            sample.primary_classification = Classifications.model_construct(
                field_uri=Classifications._convert_id_to_uri(summary["uuid"]),
                job=Jobs.model_construct(
                    field_uri=Jobs._convert_id_to_uri(summary["job_uuid"]),
                    name=summary["job_name"],
                ),
                sample=sample,
                # the endpoint only returns successful analyses
                complete=True,
                success=True,
                results_uri=summary["results_uri"],
            )

        profile = datum.get("functional_profile")
        if profile:
            functional_profiles[sample.id] = FunctionalProfiles.model_construct(
                field_uri=FunctionalProfiles._convert_id_to_uri(profile["uuid"]),
                sample=sample,
                complete=True,
                success=True,
                results_uri=profile["results_uri"],
            )

        samples.append(sample)

    return samples, functional_profiles


class SampleCollection(BaseSampleCollection):
    def __init__(
        self,
        objects: list[Samples],
        *,
        functional_profiles: dict[str, FunctionalProfiles] | None = None,
        **kwargs,
    ):
        super().__init__(objects, **kwargs)
        # Cache functional profiles for `self._functional_profiles`
        self._kwargs["functional_profiles"] = functional_profiles or {}

    @cached_property
    def _functional_profiles(self) -> list[FunctionalProfiles]:
        # The base implementation queries the API, using local cache instead
        profiles = self._kwargs["functional_profiles"]
        return [profiles[sample.id] for sample in self._res_list if sample.id in profiles]

    def plot(self, params: PlotParams) -> PlotResults:
        result = self._run_with_plot_error_handling(lambda: self._plot(params))
        if result.params is None:
            result.params = params
        return result

    def _run_with_plot_error_handling(self, plot_fn: Callable[[], PlotResults]) -> PlotResults:
        import altair as alt

        with warnings.catch_warnings(record=True) as captured_warnings:
            warnings.simplefilter("always", OneCodexUserWarning)

            try:
                result = plot_fn()
            except (ValidationError, PlottingException, NoTaxaException) as e:
                return PlotResults(error=str(e))
            except alt.MaxRowsError:
                return PlotResults(
                    error="The selected dataset is too large to plot. Please try a different plot "
                    "type or select a fewer number of samples.",
                )

        # deduplicate warning messages
        seen = set()
        for warning in captured_warnings:
            if warning.category is OneCodexUserWarning:
                message = str(warning.message)
                if message not in seen:
                    seen.add(message)
                    result.warnings.append(message)
            else:
                warnings.warn(warning.message, warning.category)

        return result

    def _plot(self, params: PlotParams) -> PlotResults:
        self._validate_plot_params(params)

        if params.filter_by and params.filter_value:
            # Create a *new* filtered SampleCollection and reassign to `self`, rather than filtering
            # the current `self` in-place. We don't want to cache the filtered SampleCollection
            self = self._filter_by_metadata(params.filter_by, params.filter_value)

        if params.metric == Metric.Auto:
            params = params.model_copy(
                update={"metric": self.automatic_metric}
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
                    metric=params.metric,
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
                    metric=params.metric,
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
                metric=params.metric,
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
                    metric=params.metric,
                    diversity_metric=params.beta_metric,
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
                    metric=params.metric,
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
                    metric=params.metric,
                    diversity_metric=params.beta_metric,
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

        return PlotResults(
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
                rows = metadata.loc[metadata["sample_id"] == sample.id, field]
                if len(rows) == 1:
                    metadatum = rows.iloc[0]
                else:
                    # ambiguous or no match
                    metadatum = None
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

    def stats(self, params: StatsParams) -> StatsResults:
        with warnings.catch_warnings():
            # Turn OneCodexUserWarning into an exception (run stats in "strict" mode)
            warnings.filterwarnings("error", category=OneCodexUserWarning)

            try:
                stats_results = self._stats(params)
            except (ValidationError, StatsException, OneCodexUserWarning, NoTaxaException) as e:
                # Expected user error
                return StatsResults(params=params, error=str(e))

        if params.stats_type == StatsType.AlphaDiversity:
            plot_params = PlotParams(
                **params.model_dump(include=set(BaseParams.model_fields)),
                plot_type=PlotType.Alpha,
                plot_repr=None,
            )
            stats_results.plot_results = self.plot(plot_params)
        elif (
            params.stats_type == StatsType.Ancombc
            and len(stats_results.ancombc_results.significant_main_results) > 0
        ):
            stats_results.plot_results = self._run_with_plot_error_handling(
                lambda: PlotResults(
                    chart=stats_results.ancombc_results.plot(return_chart=True).to_dict()
                )
            )

        return stats_results

    def _stats(self, params: StatsParams) -> StatsResults:
        self._validate_stats_params(params)

        if params.filter_by and params.filter_value:
            self = self._filter_by_metadata(params.filter_by, params.filter_value)

        if params.metric == Metric.Auto:
            params = params.model_copy(update={"metric": self.automatic_metric})

        if params.secondary_group_by:
            group_by = (params.group_by, params.secondary_group_by)
        else:
            group_by = params.group_by

        match params.stats_type:
            case StatsType.AlphaDiversity:
                alpha_diversity_results = self.alpha_diversity_stats(
                    group_by=group_by,
                    paired_by=params.paired_by,
                    metric=params.metric,
                    diversity_metric=params.alpha_metric,
                    rank=params.rank,
                )
                return StatsResults(params=params, alpha_diversity_results=alpha_diversity_results)
            case StatsType.BetaDiversity:
                beta_diversity_results = self.beta_diversity_stats(
                    group_by=group_by,
                    metric=params.metric,
                    diversity_metric=params.beta_metric,
                    rank=params.rank,
                )
                return StatsResults(params=params, beta_diversity_results=beta_diversity_results)
            case StatsType.Ancombc:
                ancombc_results = self._ancombc(
                    group_by=group_by,
                    reference_group=params.reference_group,
                    metric=params.metric,
                    rank=params.rank,
                )
                return StatsResults(params=params, ancombc_results=ancombc_results)
            case _:
                raise OneCodexException(f"Unknown stats type: {params.stats_type}")

    def _validate_stats_params(self, params: StatsParams):
        if not self._classifications:
            raise ValidationError(
                "Classification results are not available for any of the selected samples."
            )

        for attr in METADATA_FIELD_STATS_PARAMS:
            field_value = getattr(params, attr)
            if field_value is not None and field_value not in self.metadata.columns:
                attr_display_name = attr.replace("_", " ").title()
                raise ValidationError(
                    f"The metadata field {field_value!r} does not exist. Please select a valid "
                    f"metadata field in the {attr_display_name} dropdown."
                )
