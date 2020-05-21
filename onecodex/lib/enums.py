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


class AlphaDiversityMetric(BaseEnum):
    Simpson = "simpson"
    Chao1 = "chao1"
    Shannon = "shannon"


class BetaDiversityMetric(BaseEnum):
    Jaccard = "jaccard"
    BrayCurtis = "braycurtis"
    CityBlock = "cityblock"
    WeightedUnifrac = "weighted_unifrac"
    UnweightedUnifrac = "unweighted_unifrac"


class NormalizedRank(BaseEnum):
    Superkingdom = "superkingdom"
    Kingdom = "kingdom"
    Phylum = "phylum"
    Class = "class"
    Order = "order"
    Family = "family"
    Genus = "genus"
    Species = "species"


class Rank(BaseEnum):
    """Python does not allow extending enums so we have to repeat the ranks."""

    Superkingdom = "superkingdom"
    Kingdom = "kingdom"
    Phylum = "phylum"
    Class = "class"
    Order = "order"
    Family = "family"
    Genus = "genus"
    Species = "species"
    Auto = "auto"


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
