from __future__ import annotations

import io
import csv
from typing import TYPE_CHECKING

from onecodex.exceptions import OneCodexException
from .enums import PlotType, PlotRepr, ExportFormat

if TYPE_CHECKING:
    import altair as alt
    import pandas as pd

    from .models import PlotParams

METRIC_COLUMNS = [
    "Reads",
    "Mean Reads",
    "Reads (Normalized)",
    "Relative Abundance",
    "Mean Relative Abundance",
]


def export_chart_data(params: PlotParams, chart: alt.Chart) -> bytes | None:
    """Export chart data to CSV or XLSX."""
    if not params.export_format:
        return None

    df = _extract_chart_data(params, chart)
    return _dataframe_to_bytes(params, df)


def _extract_chart_data(params: PlotParams, chart: alt.Chart) -> pd.DataFrame:
    chart_data = chart.hconcat[1].data if params.plot_repr == PlotRepr.Distance else chart.data
    df = chart_data.copy()

    if "url" in df.columns:
        df.drop("url", axis=1, inplace=True)

    if params.plot_type == PlotType.Taxa:
        values_col = None
        for col_name in METRIC_COLUMNS:
            if col_name in df.columns:
                values_col = col_name

        if not values_col:
            raise OneCodexException("Cannot export chart data: unexpected or missing metric column")

        columns = "Label"
        if params.group_by:
            columns = params.group_by

        df = df.pivot_table(
            columns=columns, values=values_col, index=["tax_name"], sort=False
        ).reset_index(drop=False)
        df.rename(columns={"tax_name": "Tax Name"}, inplace=True)
    elif params.plot_type == PlotType.Alpha:
        df.rename(columns={"Label": "Sample"}, inplace=True)
        alpha_metric_col_name = " ".join([part.title() for part in params.alpha_metric.split("_")])
        df.rename(columns={params.alpha_metric: alpha_metric_col_name}, inplace=True)

        columns_to_include = []
        if params.group_by:
            columns_to_include.append(params.group_by)
        if params.secondary_group_by and params.secondary_group_by not in columns_to_include:
            columns_to_include.append(params.secondary_group_by)
        if params.facet_by and params.facet_by not in columns_to_include:
            columns_to_include.append(params.facet_by)

        if columns_to_include:
            columns_to_include = ["Sample"] + columns_to_include + [alpha_metric_col_name]
            df = df[columns_to_include]
        else:
            df.rename(columns={"classification_id": "Classification ID"}, inplace=True)
            df = df[["Sample", "Classification ID", alpha_metric_col_name]]
    elif params.plot_type == PlotType.Beta:
        if params.plot_repr == PlotRepr.Distance:
            df = (
                df.pivot(index="1) Label", columns="2) Label", values="Distance")
                .reset_index()
                .fillna(0)
                .rename(columns={"1) Label": "Sample"})
            )
        elif params.plot_repr == PlotRepr.Pca or params.plot_repr == PlotRepr.Pcoa:
            df.rename(
                columns={"Label": "Sample", "classification_id": "Classification ID"},
                inplace=True,
            )
            df = df[["Sample", "Classification ID", "PC1", "PC2"]]
    elif params.plot_type == PlotType.Functional:
        df = (
            df.pivot_table(index=["function_name"], columns=["Label"], values="value", sort=False)
            .reset_index(drop=False)
            .fillna(0)
        )
        df.rename(columns={"function_name": "Function Name"}, inplace=True)
    else:
        raise NotImplementedError(params.plot_type)

    return df


def _dataframe_to_bytes(params: PlotParams, df: pd.DataFrame) -> bytes:
    import pandas as pd

    if params.export_format == ExportFormat.Csv:
        return df.to_csv(index=False, quoting=csv.QUOTE_NONNUMERIC).encode("utf-8")
    elif params.export_format == ExportFormat.Xlsx:
        buffer = io.BytesIO()
        with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name="Custom Plots Data", index=False)
        return buffer.getvalue()
    else:
        raise NotImplementedError(params.export_format)
