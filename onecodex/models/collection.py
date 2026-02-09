from onecodex.models.base_sample_collection import BaseSampleCollection
from onecodex.stats import StatsMixin


from onecodex.viz import (
    VizBargraphMixin,
    VizDistanceMixin,
    VizFunctionalHeatmapMixin,
    VizHeatmapMixin,
    VizMetadataMixin,
    VizPCAMixin,
)


class SampleCollection(
    StatsMixin,
    VizDistanceMixin,  # < DistanceMixin < TaxonomyMixin
    VizHeatmapMixin,
    VizPCAMixin,
    VizMetadataMixin,
    VizBargraphMixin,
    VizFunctionalHeatmapMixin,
    BaseSampleCollection,
):
    """A collection of `Samples` or `Classifications` objects.

    Includes lots of methods for analysis of classifications results.
    """
