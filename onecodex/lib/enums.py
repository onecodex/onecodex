from enum import Enum

__all__ = ["Metric", "AbundanceMetric", "AlphaDiversityMetric", "BetaDiversityMetric"]


class BaseEnum(str, Enum):
    @classmethod
    def has_value(cls, value):
        return value in cls.values()

    @classmethod
    def values(cls):
        return [e.value for e in cls]


class AbundanceMetric(BaseEnum):
    Abundance = "abundance"
    AbundanceWChildren = "abundance_w_children"


class Metric(BaseEnum):
    Auto = "auto"
    Readcount = "readcount"
    ReadcountWChildren = "readcount_w_children"
    Abundance = "abundance"
    AbundanceWChildren = "abundance_w_children"

    @property
    def dtype(self):
        if self == Metric.Readcount or self == Metric.ReadcountWChildren:
            return int
        elif self == Metric.Abundance or self == Metric.AbundanceWChildren:
            return float
        else:
            raise ValueError("Metric {} does not have an associated dtype.".format(self))


class AlphaDiversityMetric(BaseEnum):
    Simpson = "simpson"
    ObservedTaxa = "observed_taxa"
    Shannon = "shannon"
    Chao1 = "chao1"


class BetaDiversityMetric(BaseEnum):
    Jaccard = "jaccard"
    BrayCurtis = "braycurtis"
    CityBlock = "cityblock"
    Manhattan = "manhattan"
    WeightedUnifrac = "weighted_unifrac"
    UnweightedUnifrac = "unweighted_unifrac"
    Aitchison = "aitchison"


class Rank(BaseEnum):
    Superkingdom = "superkingdom"
    Kingdom = "kingdom"
    Phylum = "phylum"
    Class = "class"
    Order = "order"
    Family = "family"
    Genus = "genus"
    Species = "species"
    Auto = "auto"

    @property
    def level(self):
        if self not in _RANK_TO_LEVEL:
            raise ValueError(f"Rank {self} has no level.")
        return _RANK_TO_LEVEL[self]


_RANK_TO_LEVEL = {
    Rank.Species: 0,
    Rank.Genus: 1,
    Rank.Family: 2,
    Rank.Order: 3,
    Rank.Class: 4,
    Rank.Phylum: 5,
    Rank.Kingdom: 6,
    Rank.Superkingdom: 7,
}


class Linkage(BaseEnum):
    Average = "average"
    Single = "single"
    Complete = "complete"
    Weighted = "weighted"
    Centroid = "centroid"
    Median = "median"


class OrdinationMethod(BaseEnum):
    Pcoa = "pcoa"
    Smacof = "smacof"


class AnalysisType(BaseEnum):
    Classification = "classification"
    Functional = "functional"


class FunctionalAnnotations(BaseEnum):
    Pathways = "pathways"
    MetaCyc = "metacyc"
    EggNog = "eggnog"
    Go = "go"
    Ko = "ko"
    Ec = "ec"
    Pfam = "pfam"
    Reaction = "reaction"


class FunctionalAnnotationsMetric(BaseEnum):
    Rpk = "rpk"
    Cpm = "cpm"
    Abundance = "abundance"
    Coverage = "coverage"


class Link(BaseEnum):
    Ocx = "ocx"
    Ncbi = "ncbi"
