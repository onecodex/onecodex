from onecodex.lib.enums import FunctionalAnnotations, Metric
from .models import PlotParams
from .enums import PlotType, PlotRepr
from typing import Any
import asyncio
import time
from collections import deque


class AsyncRateLimiter:
    """Simple async rate limiter: at most `max_actions` per `period` seconds."""

    def __init__(self, max_actions: int, period: float):
        self.max_actions = max_actions
        self.period = period
        self.actions = deque()
        self.lock = asyncio.Lock()

    async def acquire(self):
        async with self.lock:
            while True:
                now = time.monotonic()

                # Remove timestamps older than `period`
                cutoff = now - self.period
                while self.actions and self.actions[0] <= cutoff:
                    self.actions.popleft()

                if len(self.actions) < self.max_actions:
                    self.actions.append(now)
                    return

                # Sleep until the oldest action falls out of the window, make
                # sure sleep is not negative
                sleep_for = max(self.period - (now - self.actions[0]), 0)
                await asyncio.sleep(sleep_for)


def get_plot_title(params: PlotParams) -> str:
    source_name = params.source_name
    start = ""
    if params.facet_by:
        start = params.facet_by
    else:
        if params.plot_type == PlotType.Taxa:
            if params.metric == Metric.Readcount:
                start = "Normalized readcount"
            elif params.metric == Metric.ReadcountWChildren:
                start = "Normalized readcount with children"
            elif params.metric in {Metric.Abundance, Metric.AbundanceWChildren}:
                start = "Relative abundance"
            else:
                start = "Taxa"
        elif params.plot_type == PlotType.Alpha:
            start = "Alpha diversity"
        elif params.plot_type == PlotType.Beta:
            if params.plot_repr == PlotRepr.Distance:
                start = "Distance"
            elif params.plot_repr == PlotRepr.Pca:
                start = "PCA"
            elif params.plot_repr == PlotRepr.Pcoa:
                start = "PCoA"
        elif params.plot_type == PlotType.Functional:
            if params.functional_annotation == FunctionalAnnotations.Pathways:
                start = params.functional_pathways_metric.plot_label
            else:
                start = params.functional_metric.plot_label
            source_name = f"top {params.functional_top_n} " + source_name

    return f"{start} plot of {source_name} samples"


def format_classification_results(results: dict) -> dict[str, Any]:
    """Format the ClassificationRun results.json as a table."""

    table = []
    has_abundances = results.get("has_abundances", False)

    for tax_id, data in results.get("data", {}).items():
        entry = {
            "tax_id": str(tax_id),
            "parent_tax_id": data.get("p"),
            "readcount": data.get("n", 0),
            "readcount_w_children": data["t"],
            "name": data["i"],
            "rank": data["r"],
            "abundance": data.get("a") if has_abundances else None,
            "abundance_w_children": data.get("c") if has_abundances else None,
        }

        # Skip 0 values (i.e., reads with only k-mer hits but no reads or child reads)
        if (
            entry["readcount"] < 1
            and entry["readcount_w_children"] < 1
            and (entry["abundance"] is None or entry["abundance"] <= 0.0)
            and (entry.get("abundance_w_children") is None or entry["abundance_w_children"] <= 0.0)
        ):
            continue

        table.append(entry)

    return {
        "table": sorted(table, key=lambda r: r["readcount_w_children"], reverse=True),
        "n_reads": results["n_reads"],
        "host_tax_ids": results["hosts"],
    }
