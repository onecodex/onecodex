import math
from typing import Any, NamedTuple, Optional

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
        # no need to rehydrate rows we don't need, skip certain functional ids to maintain
        # backwards compat. with data returned from /filtered_results.
        if (
            (annotation_filter is not None and group_name != annotation_filter)
            or (metric_filter is not None and metric != metric_filter)
            or (metric_filter is not None and feature_id in _SKIP_FUNCTIONAL_IDS)
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


class _StandardContribution(NamedTuple):
    taxon_id: Any
    cpm: float
    rpk: float


class _StandardFeature(NamedTuple):
    id: str
    name: str
    cpm: float
    rpk: float
    contributions: list


class _PathwayContribution(NamedTuple):
    taxon_id: Any
    abundance: float
    coverage: float
    metacyc_cpm: float


class _PathwayFeature(NamedTuple):
    id: str
    name: str
    abundance: float
    coverage: float
    metacyc_cpm: float
    contributions: list


_STANDARD_METRIC_FIELDS = {
    FunctionalAnnotationsMetric.Cpm: "cpm",
    FunctionalAnnotationsMetric.Rpk: "rpk",
}
_PATHWAY_METRIC_FIELDS = {
    FunctionalAnnotationsMetric.Abundance: "abundance",
    FunctionalAnnotationsMetric.Coverage: "coverage",
}
_METACYC_METRIC_FIELDS = {
    # metacyc values are folded into pathways
    FunctionalAnnotationsMetric.Cpm: "metacyc_cpm",
    FunctionalAnnotationsMetric.Rpk: "abundance",
}


def _normalize_taxon_id(value: Any) -> str:
    """Coerce a taxon id to a string for use in a DataFrame index.

    Required to consistently handle missing values, ints and floats (like 386414.0).
    """

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
) -> dict:
    """Select one metric directly from condensed functional results.

    Standard functional-group rows have the form:

        [id, name, community_cpm, total_rpk, contributions]

    where each contribution is:

        [taxon_id, cpm, rpk]

    Pathway rows have the form:

        [
            id,
            name,
            community_abundance,
            community_coverage,
            metacyc_cpm,
            contributions,
        ]

    where each contribution is:

        [taxon_id, abundance, coverage, metacyc_cpm]
    """

    annotation = FunctionalAnnotations.from_value(annotation)
    metric = FunctionalAnnotationsMetric.from_value(metric)
    allowed_metrics = FunctionalAnnotationsMetric.metrics_for_annotation(annotation)

    if metric not in allowed_metrics:
        raise OneCodexException(
            f"metric {metric} cannot be retrieved for functional group {annotation}"
        )

    results_group = annotation.value
    require_complete_pathway = metric == FunctionalAnnotationsMetric.CompleteAbundance

    if annotation == FunctionalAnnotations.Pathways:
        feature_type, contribution_type = _PathwayFeature, _PathwayContribution
        value_field = _PATHWAY_METRIC_FIELDS[
            FunctionalAnnotationsMetric.Abundance if require_complete_pathway else metric
        ]
    elif annotation == FunctionalAnnotations.MetaCyc:
        results_group = FunctionalAnnotations.Pathways.value
        feature_type, contribution_type = _PathwayFeature, _PathwayContribution
        value_field = _METACYC_METRIC_FIELDS[metric]
    else:
        feature_type, contribution_type = _StandardFeature, _StandardContribution
        value_field = _STANDARD_METRIC_FIELDS[metric]

    feature_ids: list[str] = []
    values: list[float] = []
    feature_name_map: dict[str, str | None] = {}

    taxon_ids: list[str] | None = [] if taxa_stratified else None

    features = condensed_results["results"].get(results_group, [])

    for row in features:
        feature = feature_type._make(row)
        feature_id = str(feature.id)
        encoded_name = feature.name

        if feature_id in _SKIP_FUNCTIONAL_IDS:
            continue

        # complete_abundance includes only pathways whose community-level
        # coverage is exactly 1.0. The reported value is still abundance.
        if require_complete_pathway and feature.coverage != 1.0:
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
            values.append(getattr(feature, value_field))
            feature_name_map[feature_id] = feature_name
            continue

        assert taxon_ids is not None

        for raw_contribution in feature.contributions:
            contribution = contribution_type._make(raw_contribution)
            taxon_id = _normalize_taxon_id(contribution.taxon_id)

            feature_ids.append(feature_id)
            taxon_ids.append(taxon_id)
            values.append(getattr(contribution, value_field))

            # Only add names for observations that were actually emitted. This
            # keeps the feature-name map aligned with dataframe columns when a
            # feature has no taxonomic contributions.
            feature_name_map[feature_id] = feature_name

    return {
        "feature_ids": feature_ids,
        "values": values,
        "feature_name_map": feature_name_map,
        "taxon_ids": taxon_ids,
        "n_reads": condensed_results["n_reads"],
        "n_mapped": condensed_results["n_mapped"],
    }
