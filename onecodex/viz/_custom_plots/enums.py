from onecodex.lib.enums import BaseEnum


class PlotType(BaseEnum):
    Taxa = "taxa"
    Alpha = "alpha"
    Beta = "beta"
    Functional = "functional"


class PlotRepr(BaseEnum):
    Bargraph = "bargraph"
    Pcoa = "pcoa"
    Pca = "pca"
    Distance = "distance"
    Heatmap = "heatmap"


class ExportFormat(BaseEnum):
    Csv = "csv"
    Xlsx = "xlsx"


class SuggestionType(BaseEnum):
    Project = "project"
    Tag = "tag"


class SamplesFilter(BaseEnum):
    WithClassifications = "with_classifications"
    WithFunctionalResults = "with_functional_results"
