from onecodex.lib.enums import FunctionalAnnotations, Metric
from .models import PlotParams
from .enums import PlotType, PlotRepr


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
            elif params.plot_repr == PlotRepr.PCA:
                start = "PCA"
            elif params.plot_repr == PlotRepr.PCoA:
                start = "PCoA"
        elif params.plot_type == PlotType.Functional:
            if params.functional_annotation == FunctionalAnnotations.Pathways:
                start = params.functional_pathways_metric.plot_label
            else:
                start = params.functional_metric.plot_label
            source_name = f"top {params.functional_top_n} " + source_name

    return f"{start} plot of {source_name} samples"
