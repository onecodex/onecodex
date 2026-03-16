from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from pydantic import BaseModel, ConfigDict, field_validator

from onecodex.lib.enums import (
    FunctionalAnnotations,
    Metric,
    FunctionalLabel,
    FunctionalAnnotationsMetric,
    Rank,
    AlphaDiversityMetric,
    BetaDiversityMetric,
)
from onecodex.stats import AlphaDiversityStatsResults, BetaDiversityStatsResults, PosthocResults

if TYPE_CHECKING:
    import pandas as pd

from .enums import ExportFormat, PlotType, PlotRepr, StatsType


class BaseParams(BaseModel):
    model_config = ConfigDict(extra="forbid", arbitrary_types_allowed=False)

    tag: str | None = None
    project: str | None = None

    metric: Metric = Metric.Auto
    alpha_metric: AlphaDiversityMetric = AlphaDiversityMetric.Shannon
    beta_metric: BetaDiversityMetric = BetaDiversityMetric.BrayCurtis
    rank: Rank = Rank.Species

    secondary_group_by: str | None = None
    filter_by: str | None = None
    filter_value: list[str] = []


class PlotParams(BaseParams):
    source_name: str = ""
    plot_type: PlotType = PlotType.Taxa
    plot_repr: PlotRepr | None = PlotRepr.Bargraph
    export_format: ExportFormat | None = None

    top_n: int = 10
    facet_by: str | None = None
    group_by: str | None = None
    label_by: list[str] = []
    sort_by: str | None = None

    functional_annotation: FunctionalAnnotations = FunctionalAnnotations.Go
    functional_metric: FunctionalAnnotationsMetric = FunctionalAnnotationsMetric.Cpm
    functional_pathways_metric: FunctionalAnnotationsMetric = FunctionalAnnotationsMetric.Abundance
    functional_top_n: int = 10
    functional_label: FunctionalLabel = FunctionalLabel.Name

    @field_validator("label_by")
    def validate_label_by(cls, value):
        if not value:
            return ["sample_name"]
        return value


@dataclass(kw_only=True)
class PlotResult:
    params: PlotParams
    chart: dict | None = None
    x_axis_label_links: dict[str, str] = field(default_factory=dict)
    error: str | None = None
    warnings: list[str] = field(default_factory=list)
    exported_chart_data: bytes | None = None

    def to_dict(self) -> dict:
        return {
            "params": self.params.model_dump(),
            "chart": self.chart,
            "x_axis_label_links": self.x_axis_label_links,
            "error": self.error,
            "warnings": self.warnings,
            "exported_chart_data": self.exported_chart_data,
        }


class StatsParams(BaseParams):
    stats_type: StatsType = StatsType.AlphaDiversity
    group_by: str
    paired_by: str | None = None


@dataclass(kw_only=True)
class StatsResult:
    params: StatsParams
    results: AlphaDiversityStatsResults | BetaDiversityStatsResults | None = None
    error: str | None = None

    def to_dict(self) -> dict:
        return {
            "params": self.params.model_dump(),
            "results": _stats_results_to_dict(self.results) if self.results is not None else None,
            "error": self.error,
        }


def _stats_results_to_dict(
    results: AlphaDiversityStatsResults | BetaDiversityStatsResults,
) -> dict:
    result = {
        "test": str(results.test),
        "statistic": results.statistic,
        "pvalue": results.pvalue,
        "alpha": results.alpha,
        "sample_size": results.sample_size,
        "group_by_variable": results.group_by_variable,
        "group_sizes": results.group_sizes,
    }

    if isinstance(results, AlphaDiversityStatsResults):
        result["paired_by_variable"] = results.paired_by_variable
    elif isinstance(results, BetaDiversityStatsResults):
        result["num_permutations"] = results.num_permutations

    if results.posthoc is not None:
        result["posthoc"] = _posthoc_to_dict(results.posthoc)

    return result


def _posthoc_to_dict(posthoc: PosthocResults) -> dict:
    result = {
        "test": str(posthoc.test),
        "adjustment_method": str(posthoc.adjustment_method),
        "adjusted_pvalues": _posthoc_df_to_dict(posthoc.adjusted_pvalues),
    }
    if posthoc.pvalues is not None:
        result["pvalues"] = _posthoc_df_to_dict(posthoc.pvalues)
    if posthoc.statistics is not None:
        result["statistics"] = _posthoc_df_to_dict(posthoc.statistics)
    return result


def _posthoc_df_to_dict(df: pd.DataFrame) -> dict[str, dict[str, float]]:
    """Convert a posthoc pandas DataFrame to a nested dict mapping group pairs to values."""
    return {
        str(row_name): {str(col_name): value for col_name, value in row.items()}
        for row_name, row in df.iterrows()
    }
