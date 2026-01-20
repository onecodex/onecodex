from onecodex.exceptions import OneCodexException

try:
    # Python 3.11 changed `__str__`/`__format__` behavior for enums with mixed-in data types
    from enum import ReprEnum as Enum
except ImportError:
    from enum import Enum

__all__ = ["Metric", "AlphaDiversityMetric", "BetaDiversityMetric"]

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


class Metric(BaseEnum):
    r"""Metrics for taxonomic abundance data.

    .. attribute:: Readcount

       Number of reads assigned to a given taxon.

    .. attribute:: ReadcountWChildren

       Read count for a taxon and all its taxonomic descendants.

    .. attribute:: PropReadcount

       Readcount as a proportion of total reads in the sample.

       .. math::

          \frac{\text{readcount}}{\text{n\_reads\_total}}

    .. attribute:: PropReadcountWChildren

       ReadcountWChildren as a proportion of total reads in the sample.

    .. attribute:: PropClassified

       Readcount as a proportion of classified microbial reads.

       .. math::

          \frac{\text{readcount}}{\text{n\_mapped\_microbial\_reads}}

       Where ``n_mapped_microbial_reads`` excludes host reads (typically human/mouse)
       and nonspecific reads (root and cellular organisms).

    .. attribute:: PropClassifiedWChildren

       ReadcountWChildren as a proportion of classified microbial reads.

    .. attribute:: Abundance

       Relative abundance estimate for a given taxon, computed by the classifier.

    .. attribute:: AbundanceWChildren

       Abundance for a taxon and all its descendants.

    .. attribute:: NormalizedReadcount

       Readcount normalized by the sum of readcounts across all taxa in the sample.

       .. math::

          \frac{\text{readcount}}{\sum_{\text{taxa}} \text{readcount}}

    .. attribute:: NormalizedReadcountWChildren

       NormalizedReadcount for a taxon and all its descendants.
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
    NormalizedReadcount = "normalized_readcount"
    NormalizedReadcountWChildren = "normalized_readcount_w_children"

    @property
    def is_abundance_metric(self) -> bool:
        """Returns True if this metric is an abundance metric (Abundance or AbundanceWChildren)."""
        return self in {Metric.Abundance, Metric.AbundanceWChildren}

    @property
    def is_normalized(self) -> bool:
        """Return true if the metric has been normalized (ie proportionalized) in some way."""
        return {
            Metric.Abundance: True,
            Metric.AbundanceWChildren: True,
            Metric.PropReadcount: True,
            Metric.PropReadcountWChildren: True,
            Metric.PropClassified: True,
            Metric.PropClassifiedWChildren: True,
            Metric.NormalizedReadcount: True,
            Metric.NormalizedReadcountWChildren: True,
            Metric.Readcount: False,
            Metric.ReadcountWChildren: False,
        }[self]

    @property
    def includes_children(self) -> bool:
        return self in (
            Metric.ReadcountWChildren,
            Metric.AbundanceWChildren,
            Metric.PropReadcountWChildren,
            Metric.NormalizedReadcountWChildren,
            Metric.PropClassifiedWChildren,
        )

    @property
    def results_key(self):
        """Return the key used to fetch the raw value for this metric in Classifications.results."""
        return {
            Metric.Readcount: "readcount",
            Metric.ReadcountWChildren: "readcount_w_children",
            Metric.PropReadcount: "readcount",
            Metric.PropReadcountWChildren: "readcount_w_children",
            Metric.PropClassified: "readcount",
            Metric.PropClassifiedWChildren: "readcount_w_children",
            Metric.Abundance: "abundance",
            Metric.AbundanceWChildren: "abundance_w_children",
            Metric.NormalizedReadcount: "readcount",
            Metric.NormalizedReadcountWChildren: "readcount_w_children",
        }[self]

    @property
    def dtype(self):
        return {
            Metric.Readcount: int,
            Metric.ReadcountWChildren: int,
            Metric.Abundance: float,
            Metric.AbundanceWChildren: float,
            Metric.PropReadcount: float,
            Metric.PropReadcountWChildren: float,
            Metric.PropClassified: float,
            Metric.PropClassifiedWChildren: float,
            Metric.NormalizedReadcount: float,
            Metric.NormalizedReadcountWChildren: float,
        }[self]

    @property
    def display_name(self) -> str:
        return {
            Metric.Readcount: "Readcount",
            Metric.ReadcountWChildren: "Readcount With Children",
            Metric.Abundance: "Relative Abundance",
            Metric.AbundanceWChildren: "Relative Abundance",
            Metric.PropReadcount: "Reads (Normalized)",
            Metric.PropReadcountWChildren: "Reads (Normalized)",
            Metric.PropClassified: "Classified Reads (Normalized)",
            Metric.PropClassifiedWChildren: "Classified Reads (Normalized)",
            Metric.NormalizedReadcount: "Normalized Readcount",
            Metric.NormalizedReadcountWChildren: "Normalized Readcount With Children",
        }[self]


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
