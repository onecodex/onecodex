import math
from dataclasses import dataclass
from typing import Any, Optional

from onecodex.exceptions import OneCodexException
from onecodex.lib.enums import FunctionalAnnotations, FunctionalAnnotationsMetric

_SKIP_FUNCTIONAL_IDS = {"UNMAPPED", "UNGROUPED", "UNINTEGRATED"}


def _rehydrate_functional_results(
    condensed_results: dict,
    annotation_filter: Optional[str] = None,
    metric_filter: Optional[str] = None,
    taxa_stratified_filter: Optional[bool] = None,
) -> dict:
    """Rehydrate condensed functional results into the public API results format.

    the annotation, metric, and taxa stratified filters can be used to control what gets
    rehydrated. If any of these are set, all other annotations/metrics/stratifications will
    be removed prior to returning the rehydrated results.
    """

    # maps tax IDs -> names
    taxa_map = {node["id"]: node.get("name") for node in condensed_results["taxonomy"]["nodes"]}
    taxa_map["0"] = "unclassified"

    table = []
    pathway_table = []

    # do we only want things with complete (1.0) abundance?
    complete_abundance_filter = metric_filter == FunctionalAnnotationsMetric.CompleteAbundance

    if complete_abundance_filter:
        metric_filter = FunctionalAnnotationsMetric.Abundance.value

    def add_row(
        group_name: str,
        feature_id: str,
        feature_name: Optional[str],
        metric: str,
        value: float,
        taxa_stratified: bool = False,
        taxon_id: Optional[str] = None,
        destination: Optional[list] = None,
    ) -> None:
        """Add a row formatted for results['table']."""
        # no need to rehydrate rows we don't need
        if (
            (annotation_filter is not None and group_name != annotation_filter)
            or (metric_filter is not None and metric != metric_filter)
            or (taxa_stratified_filter is not None and taxa_stratified != taxa_stratified_filter)
        ):
            return

        destination = table if destination is None else destination
        destination.append(
            {
                "group_name": group_name,
                "id": feature_id,
                "name": feature_name,
                "metric": metric,
                "value": value,
                "taxa_stratified": taxa_stratified,
                "taxon_id": taxon_id,
                "taxon_name": taxa_map.get(taxon_id) if taxa_stratified else None,
            }
        )

    def add_species_rows(
        group_name: str,
        feature_id: str,
        feature_name: Optional[str],
        taxon_id: str,
        metric_values: tuple[tuple[str, float], ...],
        destination: Optional[list] = None,
    ) -> None:
        """Add a species-level (taxa stratified) row formatted for results['table'].

        For all functional groups, except pathways, a taxa stratified row will have one of
        two metrics: cpm or rpk. For pathways, these metrics are abundance (synonymous with RPK)
        and coverage.
        """
        for metric, value in metric_values:
            add_row(
                group_name,
                feature_id,
                feature_name,
                metric,
                value,
                taxa_stratified=True,
                taxon_id=taxon_id,
                destination=destination,
            )

    # standard functional groups use the following format:
    # [id, name, total_cpm, total_rpk, [[taxid, cpm, rpk], ...]]
    for group_name, features in condensed_results["results"].items():
        # we'll do pathways separately because we neeed to split pathways into
        # the metacyc functional group, has cpm/rpk metrics, and the pathways
        # group, which has abundance/coverage metrics
        if group_name == "pathways":
            continue

        for feature_id, encoded_name, community_cpm, total_rpk, contributions in features:
            # Missing metadata is encoded by repeating the feature ID.
            feature_name = None if encoded_name == feature_id else encoded_name

            if taxa_stratified_filter is not True:
                add_row(group_name, feature_id, feature_name, "cpm", community_cpm)
                add_row(group_name, feature_id, feature_name, "rpk", total_rpk)

            if taxa_stratified_filter is not False:
                for taxon_id, cpm, rpk in contributions:
                    add_species_rows(
                        group_name,
                        feature_id,
                        feature_name,
                        taxon_id,
                        (("rpk", rpk), ("cpm", cpm)),
                    )

    pathways = condensed_results["results"].get("pathways", [])

    # metacyc values are folded into each condensed pathway:
    # [id, name, community_abundance, community_coverage, metacyc_cpm,
    #  [[taxid, species_abundance, species_coverage, metacyc_cpm], ...]]
    for (
        pathway_id,
        pathway_name,
        community_abundance,
        community_coverage,
        community_cpm,
        contributions,
    ) in pathways:
        if complete_abundance_filter and (community_coverage != 1.0):
            continue

        # these are only added when not stratifying by taxa
        if taxa_stratified_filter is not True:
            add_row("metacyc", pathway_id, pathway_name, "cpm", community_cpm)
            add_row("metacyc", pathway_id, pathway_name, "rpk", community_abundance)

            add_row(
                "pathways",
                pathway_id,
                pathway_name,
                "coverage",
                community_coverage,
                destination=pathway_table,
            )
            add_row(
                "pathways",
                pathway_id,
                pathway_name,
                "abundance",
                community_abundance,
                destination=pathway_table,
            )

        if taxa_stratified_filter is not False:
            for taxon_id, abundance, coverage, cpm in contributions:
                # metacyc ids/names are synonymous with pathway ids/names
                add_species_rows(
                    "metacyc",
                    pathway_id,
                    pathway_name,
                    taxon_id,
                    (("rpk", abundance), ("cpm", cpm)),
                )
                # add these to pathways separately so they aren't interleaved with metacyc
                add_species_rows(
                    "pathways",
                    pathway_id,
                    pathway_name,
                    taxon_id,
                    (("coverage", coverage), ("abundance", abundance)),
                    destination=pathway_table,
                )

    table.extend(pathway_table)

    return {
        "table": table,
        "n_reads": condensed_results["n_reads"],
        "n_mapped": condensed_results["n_mapped"],
    }


_SKIP_FUNCTIONAL_IDS = {"UNMAPPED", "UNGROUPED", "UNINTEGRATED"}


@dataclass(frozen=True)
class FunctionalResultValues:
    """A selected metric from one functional profile.

    ``feature_ids``, ``values``, and—when stratified—``taxon_ids`` are
    positionally aligned.

    Names are stored as maps because they repeat across many observations.
    """

    feature_ids: list[str]
    values: list[float]
    feature_name_map: dict[str, str | None]
    taxon_ids: list[str] | None
    taxon_name_map: dict[str, str | None]
    n_reads: int
    n_mapped: int

    def __post_init__(self) -> None:
        if len(self.feature_ids) != len(self.values):
            raise ValueError("feature_ids and values must have equal lengths")

        if self.taxon_ids is not None and len(self.taxon_ids) != len(self.values):
            raise ValueError("taxon_ids and values must have equal lengths")


@dataclass(frozen=True)
class _SelectionSpec:
    """Positions for one metric in the condensed representation."""

    results_group: str
    community_value_index: int
    contribution_value_index: int
    require_complete_pathway: bool = False


_STANDARD_METRIC_INDEXES = {
    FunctionalAnnotationsMetric.Cpm: (2, 1),
    FunctionalAnnotationsMetric.Rpk: (3, 2),
}

_PATHWAY_METRIC_INDEXES = {
    FunctionalAnnotationsMetric.Abundance: (2, 1),
    FunctionalAnnotationsMetric.Coverage: (3, 2),
}

_METACYC_METRIC_INDEXES = {
    # MetaCyc values are folded into the condensed pathway rows.
    FunctionalAnnotationsMetric.Cpm: (4, 3),
    FunctionalAnnotationsMetric.Rpk: (2, 1),
}


def _selection_spec(
    annotation: FunctionalAnnotations,
    metric: FunctionalAnnotationsMetric,
) -> _SelectionSpec:
    """Describe where a requested metric lives in the condensed arrays."""

    allowed_metrics = FunctionalAnnotationsMetric.metrics_for_annotation(annotation)
    if metric not in allowed_metrics:
        raise OneCodexException(
            f"metric {metric} cannot be retrieved for functional group {annotation}"
        )

    if annotation == FunctionalAnnotations.Pathways:
        if metric == FunctionalAnnotationsMetric.CompleteAbundance:
            community_index, contribution_index = _PATHWAY_METRIC_INDEXES[
                FunctionalAnnotationsMetric.Abundance
            ]
            return _SelectionSpec(
                results_group=FunctionalAnnotations.Pathways.value,
                community_value_index=community_index,
                contribution_value_index=contribution_index,
                require_complete_pathway=True,
            )

        community_index, contribution_index = _PATHWAY_METRIC_INDEXES[metric]
        return _SelectionSpec(
            results_group=FunctionalAnnotations.Pathways.value,
            community_value_index=community_index,
            contribution_value_index=contribution_index,
        )

    if annotation == FunctionalAnnotations.MetaCyc:
        community_index, contribution_index = _METACYC_METRIC_INDEXES[metric]
        return _SelectionSpec(
            results_group=FunctionalAnnotations.Pathways.value,
            community_value_index=community_index,
            contribution_value_index=contribution_index,
        )

    community_index, contribution_index = _STANDARD_METRIC_INDEXES[metric]
    return _SelectionSpec(
        results_group=annotation.value,
        community_value_index=community_index,
        contribution_value_index=contribution_index,
    )


def _normalize_taxon_id(value: Any) -> str:
    """Normalize taxon IDs consistently for dataframe column keys."""

    if value is None:
        return ""

    if isinstance(value, float):
        if math.isnan(value):
            return ""
        if value.is_integer():
            return str(int(value))

    return str(value)


def _select_condensed_functional_results(
    condensed_results: dict,
    annotation: FunctionalAnnotations | str,
    metric: FunctionalAnnotationsMetric | str,
    taxa_stratified: bool,
) -> FunctionalResultValues:
    """Select one metric directly from condensed functional results.

    This does not construct legacy result-row dictionaries or an intermediate
    pandas DataFrame.

    Standard functional-group rows have the form::

        [id, name, community_cpm, total_rpk, contributions]

    where each contribution is::

        [taxon_id, cpm, rpk]

    Pathway rows have the form::

        [
            id,
            name,
            community_abundance,
            community_coverage,
            metacyc_cpm,
            contributions,
        ]

    where each contribution is::

        [taxon_id, abundance, coverage, metacyc_cpm]
    """

    annotation = FunctionalAnnotations.from_value(annotation)
    metric = FunctionalAnnotationsMetric.from_value(metric)
    spec = _selection_spec(annotation, metric)

    taxonomy_names = {
        _normalize_taxon_id(node["id"]): node.get("name")
        for node in condensed_results["taxonomy"]["nodes"]
    }
    taxonomy_names["0"] = "unclassified"

    feature_ids: list[str] = []
    values: list[float] = []
    feature_name_map: dict[str, str | None] = {}

    taxon_ids: list[str] | None = [] if taxa_stratified else None
    taxon_name_map: dict[str, str | None] = {}

    features = condensed_results["results"].get(spec.results_group, [])

    for feature in features:
        feature_id = str(feature[0])
        encoded_name = feature[1]

        if feature_id in _SKIP_FUNCTIONAL_IDS:
            continue

        # complete_abundance includes only pathways whose community-level
        # coverage is exactly 1.0. The reported value is still abundance.
        if spec.require_complete_pathway and feature[3] != 1.0:
            continue

        # Missing names for standard groups are encoded by repeating the ID.
        # Pathway names are also used for the derived MetaCyc group.
        if annotation in (
            FunctionalAnnotations.Pathways,
            FunctionalAnnotations.MetaCyc,
        ):
            feature_name = encoded_name
        else:
            feature_name = None if encoded_name == feature_id else encoded_name

        if not taxa_stratified:
            feature_ids.append(feature_id)
            values.append(feature[spec.community_value_index])
            feature_name_map[feature_id] = feature_name
            continue

        assert taxon_ids is not None

        # Contributions are always the final element of a condensed feature.
        for contribution in feature[-1]:
            taxon_id = _normalize_taxon_id(contribution[0])

            feature_ids.append(feature_id)
            taxon_ids.append(taxon_id)
            values.append(contribution[spec.contribution_value_index])

            # Only add names for observations that were actually emitted. This
            # keeps the feature-name map aligned with dataframe columns when a
            # feature has no taxonomic contributions.
            feature_name_map[feature_id] = feature_name
            taxon_name_map[taxon_id] = taxonomy_names.get(taxon_id)

    return FunctionalResultValues(
        feature_ids=feature_ids,
        values=values,
        feature_name_map=feature_name_map,
        taxon_ids=taxon_ids,
        taxon_name_map=taxon_name_map,
        n_reads=condensed_results["n_reads"],
        n_mapped=condensed_results["n_mapped"],
    )


# def _functional_values(
#    self,
#    annotation: FunctionalAnnotations | str,
#    metric: FunctionalAnnotationsMetric | str,
#    taxa_stratified: bool,
# ) -> FunctionalResultValues:
#    return _select_condensed_functional_results(
#        self._condensed_results(),
#        annotation=annotation,
#        metric=metric,
#        taxa_stratified=taxa_stratified,
#    )
