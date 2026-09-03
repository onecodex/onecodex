from typing import Optional

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


def _rehydrate_filtered_functional_results(
    condensed_results: dict,
    annotation: FunctionalAnnotations,
    metric: FunctionalAnnotationsMetric,
    taxa_stratified: bool,
) -> dict:
    """Rehydrate condensed functional results into the public filtered API results format."""

    _SKIP_FUNCTIONAL_IDS = {"UNMAPPED", "UNGROUPED", "UNINTEGRATED"}

    annotation = FunctionalAnnotations.from_value(annotation)
    metric = FunctionalAnnotationsMetric.from_value(metric)

    if metric not in FunctionalAnnotationsMetric.metrics_for_annotation(annotation):
        raise OneCodexException(
            f"metric {metric} cannot be retrieved for functional group {annotation}"
        )

    effective_metric = metric

    if annotation in (FunctionalAnnotations.Pathways, FunctionalAnnotations.MetaCyc):
        selected = [
            pathway
            for pathway in condensed_results["results"].get("pathways", [])
            if pathway[0] not in _SKIP_FUNCTIONAL_IDS
        ]

        # complete abundance means annotations where the community level coverage is 1.0
        if (
            annotation == FunctionalAnnotations.Pathways
            and metric == FunctionalAnnotationsMetric.CompleteAbundance
        ):
            # pathway[3] == community_coverage
            selected = [pathway for pathway in selected if pathway[3] == 1.0]
            effective_metric = FunctionalAnnotationsMetric.Abundance

        selected_results = {"pathways": selected}
    else:
        selected_results = {
            annotation.value: [
                feature
                for feature in condensed_results["results"].get(annotation.value, [])
                if feature[0] not in _SKIP_FUNCTIONAL_IDS
            ]
        }

    filtered_condensed_results = {
        **condensed_results,
        "results": selected_results,
    }

    hydrated = _rehydrate_functional_results(
        filtered_condensed_results,
        annotation_filter=annotation.value,
        metric_filter=effective_metric.value,
        taxa_stratified_filter=taxa_stratified,
    )

    table = []

    for row in hydrated["table"]:
        filtered_row = {
            "id": row["id"],
            "name": row["name"],
            "value": row["value"],
        }

        if taxa_stratified:
            filtered_row["taxon_id"] = row["taxon_id"]
            filtered_row["taxon_name"] = row["taxon_name"]

        table.append(filtered_row)

    return {
        "table": table,
        "n_reads": condensed_results["n_reads"],
        "n_mapped": condensed_results["n_mapped"],
    }
