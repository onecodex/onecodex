from __future__ import annotations

from dataclasses import dataclass, field

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
from .enums import ExportFormat, PlotType, PlotRepr


class PlotParams(BaseModel):
    model_config = ConfigDict(extra="forbid", arbitrary_types_allowed=False)

    tag: str | None
    project: str | None
    source_name: str

    plot_type: PlotType
    plot_repr: PlotRepr | None
    metric: Metric
    alpha_metric: AlphaDiversityMetric
    beta_metric: BetaDiversityMetric
    export_format: ExportFormat | None

    top_n: int
    rank: Rank
    facet_by: str | None
    group_by: str | None
    secondary_group_by: str | None
    filter_by: str | None
    filter_value: list[str]
    label_by: list[str]
    sort_by: str | None

    functional_annotation: FunctionalAnnotations
    functional_metric: FunctionalAnnotationsMetric
    functional_pathways_metric: FunctionalAnnotationsMetric
    functional_top_n: int
    functional_label: FunctionalLabel

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
            "params": self.params.dict(),
            "chart": self.chart,
            "x_axis_label_links": self.x_axis_label_links,
            "error": self.error,
            "warnings": self.warnings,
            "exported_chart_data": self.exported_chart_data,
        }
