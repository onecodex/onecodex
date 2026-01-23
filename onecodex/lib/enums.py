from onecodex.exceptions import OneCodexException

from enum import Enum

__all__ = ["Metric", "AlphaDiversityMetric", "BetaDiversityMetric"]

from typing import TypeVar


T = TypeVar("T", bound="BaseEnum")


class BaseEnum(str, Enum):
    def __repr__(self) -> str:
        return str(self.value)

    def __str__(self) -> str:
        return str(self.value)

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

    Taxonomic descendants refer to all taxa below a given taxon in the taxonomic hierarchy
    (superkingdom > phylum > class > order > family > genus > species > strain). For example,
    if a read maps to the species Escherichia coli, that read contributes to the readcount
    for E. coli as well as the ``ReadcountWChildren`` for its genus (Escherichia), family
    (Enterobacteriaceae), and all ancestor taxa up the hierarchy.

    .. attribute:: Auto

       Determine appropriate metric automatically (see :meth:`~onecodex.models.collection.SampleCollection.to_classification_df`).

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

       Where ``n_mapped_microbial_reads = n_mapped_reads - n_host_reads - n_nonspecific_reads``.
       Host reads are those mapping to detected host organisms (typically tax IDs 9606 for human
       and 10090 for mouse). Nonspecific reads are those mapping to tax IDs 1 (root) and 131567
       (cellular organisms).

    .. attribute:: PropClassifiedWChildren

       ReadcountWChildren as a proportion of classified microbial reads.

    .. attribute:: Abundance

       Relative abundance estimate for a given taxon, computed by the classifier.

    .. attribute:: AbundanceWChildren

       Abundance for a taxon and all its descendants.

    .. attribute:: NormalizedReadcount

       Readcount normalized by the sum of readcounts for taxa at the specified rank.
       Values sum to 1.0 across taxa at that rank.

       .. math::

          \frac{\text{readcount}}{\sum_{\text{taxa at rank}} \text{readcount}}

    .. attribute:: NormalizedReadcountWChildren

       ReadcountWChildren normalized by the sum of ReadcountWChildren for taxa at the
       specified rank. Represents the proportion of reads that classified to the specified
       rank or below. Values sum to 1.0 across taxa at that rank.

       .. math::

          \frac{\text{readcount\_w\_children}}{\sum_{\text{taxa at rank}} \text{readcount\_w\_children}}
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
        return self in {
            Metric.Abundance,
            Metric.AbundanceWChildren,
            Metric.PropReadcount,
            Metric.PropReadcountWChildren,
            Metric.PropClassified,
            Metric.PropClassifiedWChildren,
            Metric.NormalizedReadcount,
            Metric.NormalizedReadcountWChildren,
        }

    @property
    def includes_children(self) -> bool:
        """Returns True if this metric aggregates over its own taxonomic descendants (children)."""
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
        """Returns the Python data type used to store this metric in a Numpy array or Pandas DataFrame."""
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
        """Returns the human-readable name used for this metric. Useful for plotting and exports."""
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
    """Supported alpha-diversity metrics."""

    Chao1 = "chao1"
    ObservedTaxa = "observed_taxa"
    Shannon = "shannon"
    Simpson = "simpson"


class BetaDiversityMetric(BaseEnum):
    """Supported beta-diversity metrics."""

    Aitchison = "aitchison"
    BrayCurtis = "braycurtis"
    CityBlock = "cityblock"
    Euclidean = "euclidean"
    Jaccard = "jaccard"
    Manhattan = "manhattan"
    UnweightedUnifrac = "unweighted_unifrac"
    WeightedUnifrac = "weighted_unifrac"


class Rank(BaseEnum):
    """A taxonomic rank."""

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
    """Clustering linkages used for distance plots."""

    Average = "average"
    Single = "single"
    Complete = "complete"
    Weighted = "weighted"
    Centroid = "centroid"
    Median = "median"


class OrdinationMethod(BaseEnum):
    """Ordination methods supported by distance plotting methods."""

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
    """Analysis types supported by :meth:`~onecodex.models.collection.SampleCollection.to_df`."""

    Classification = "classification"
    Functional = "functional"


class FunctionalAnnotations(BaseEnum):
    """Types of functional annotations."""

    Pathways = "pathways"
    MetaCyc = "metacyc"
    EggNog = "eggnog"
    Go = "go"
    Ko = "ko"
    Ec = "ec"
    Pfam = "pfam"
    Reaction = "reaction"


class FunctionalAnnotationsMetric(BaseEnum):
    """Metrics used to quantify functional annotations."""

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
