from onecodex.lib.enums import BaseEnum


class PlotType(BaseEnum):
    Taxa = "taxa"
    Alpha = "alpha"
    Beta = "beta"
    Functional = "functional"


class PlotRepr(BaseEnum):
    Bargraph = "bargraph"
    PCoA = "pcoa"
    PCA = "pca"
    Distance = "distance"
    Heatmap = "heatmap"


class ExportFormat(BaseEnum):
    Csv = "csv"
    Xlsx = "xlsx"


class SuggestionType(BaseEnum):
    Project = "Project"
    Tag = "Tag"
