from onecodex.exceptions import OneCodexException

try:
    # Python 3.11 changed `__str__`/`__format__` behavior for enums with mixed-in data types
    from enum import ReprEnum as Enum
except ImportError:
    from enum import Enum

__all__ = ["Metric", "AbundanceMetric", "AlphaDiversityMetric", "BetaDiversityMetric"]

from typing import TypeVar


T = TypeVar("T", bound="BaseEnum")


class BaseEnum(str, Enum):
    @classmethod
    def has_value(cls, value):
        return value in cls.values()

    @classmethod
    def values(cls):
        return [e.value for e in cls]

    @classmethod
    def from_value(cls: type[T], val: str) -> T:
        try:
            return cls(val)
        except ValueError:
            raise OneCodexException(f"{val} is not valid value for {cls.__name__}")


class AbundanceMetric(BaseEnum):
    Abundance = "abundance"
    AbundanceWChildren = "abundance_w_children"


class Metric(BaseEnum):
    """Metrics for taxonomic abundance data.

    Readcount: Number of reads assigned to a given taxon.
    ReadcountWChildren: Read count for a taxon and all its taxonomic descendants.
    PropReadcount: Readcount normalized by the number of classified microbial reads.
    PropReadcountWChildren: ReadcountWChildren normalized by the number of classified microbial reads.
    Abundance: Relative abundance estimate for a given taxon.
    AbundanceWChildren: Abundance for a taxon and all its descendants.
    """

    Auto = "auto"
    Readcount = "readcount"
    ReadcountWChildren = "readcount_w_children"
    PropClassified = "prop_classified"
    PropClassifiedWChildren = "prop_classified_w_children"
    PropReadcount = "prop_readcount"
    PropReadcountWChildren = "prop_readcount_w_children"
    Abundance = "abundance"
    AbundanceWChildren = "abundance_w_children"

    @property
    def dtype(self):
        dtype = {
            "readcount": int,
            "readcount_w_children": int,
            "abundance": float,
            "abundance_w_children": float,
            "prop_readcount": float,
            "prop_readcount_w_children": float,
            "prop_classified": float,
            "prop_classified_w_children": float,
        }.get(self.value)
        if dtype is None:
            raise ValueError(f"Metric {self} does not have an associated dtype.")
        return dtype

    @property
    def display_name(self) -> str:
        return {
            "readcount": "Reads",
            "readcount_w_children": "Reads",
            "abundance": "Relative Abundance",
            "abundance_w_children": "Relative Abundance",
            "prop_readcount": "Reads (Normalized)",
            "prop_readcount_w_children": "Reads (Normalized)",
            "prop_classified": "Classified Reads (Normalized)",
            "prop_classified_w_children": "Classified Reads (Normalized)",
        }[self.value]


# metrics that have been "normalized" (i.e., they are proportional)
NORMALIZED_METRICS = {
    Metric.Abundance,
    Metric.AbundanceWChildren,
    Metric.PropReadcount,
    Metric.PropReadcountWChildren,
    Metric.PropClassified,
    Metric.PropClassifiedWChildren,
}


class AlphaDiversityMetric(BaseEnum):
    Chao1 = "chao1"
    ObservedTaxa = "observed_taxa"
    Shannon = "shannon"
    Simpson = "simpson"


class BetaDiversityMetric(BaseEnum):
    Aitchison = "aitchison"
    BrayCurtis = "braycurtis"
    CityBlock = "cityblock"
    Euclidean = "euclidean"
    Jaccard = "jaccard"
    Manhattan = "manhattan"
    UnweightedUnifrac = "unweighted_unifrac"
    WeightedUnifrac = "weighted_unifrac"


class Rank(BaseEnum):
    Auto = "auto"
    Class = "class"
    Family = "family"
    Genus = "genus"
    Kingdom = "kingdom"
    Order = "order"
    Phylum = "phylum"
    Species = "species"
    Superkingdom = "superkingdom"

    @property
    def level(self):
        if self not in _RANK_TO_LEVEL:
            raise ValueError(f"Rank {self} has no level.")
        return _RANK_TO_LEVEL[self]

    @classmethod
    def from_value(cls: type[T], val: str) -> T:
        for member in cls:
            if member.value == val:
                return member
        raise OneCodexException(f"{val} is not valid value for {cls.__name__}")


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


class AlphaDiversityStatsTest(BaseEnum):
    Auto = "auto"
    Wilcoxon = "wilcoxon"
    Mannwhitneyu = "mannwhitneyu"
    Kruskal = "kruskal"


class BetaDiversityStatsTest(BaseEnum):
    Permanova = "permanova"


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
    CompleteAbundance = "complete_abundance"
    Coverage = "coverage"

    @classmethod
    def metrics_for_annotation(cls, annotation):
        return (
            [
                FunctionalAnnotationsMetric.Abundance,
                FunctionalAnnotationsMetric.CompleteAbundance,
                FunctionalAnnotationsMetric.Coverage,
            ]
            if annotation == FunctionalAnnotations.Pathways
            else [FunctionalAnnotationsMetric.Cpm, FunctionalAnnotationsMetric.Rpk]
        )

    @property
    def plot_label(self):
        if self in [FunctionalAnnotationsMetric.Cpm, FunctionalAnnotationsMetric.Rpk]:
            return self.value.upper()

        parts = self.value.split("_")
        return " ".join(x.title() for x in parts)


class FunctionalLabel(BaseEnum):
    Name = "name"
    Id = "id"


class Link(BaseEnum):
    Ocx = "ocx"
    Ncbi = "ncbi"
