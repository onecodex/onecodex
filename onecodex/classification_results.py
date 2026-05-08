"""Analysis utilities for One Codex classification results."""

from __future__ import annotations

import json
import math
from typing import Generator, Iterable, NamedTuple, TypedDict


class TaxonomyNode(NamedTuple):
    tax_id: str
    name: str
    rank: str
    parent_tax_id: str | None


# Increment this whenever the output format of summarize_analysis_json changes.
# Consumers who depend on a stable output format can use this to detect breaking changes.
RESULTS_FORMAT_VERSION = 9


TaxonEntry = TypedDict(
    "TaxonEntry",
    {
        "i": str,  # taxon name
        "r": str,  # taxon rank
        "n": int,  # readcount (absent when 0)
        "t": int,  # cumulative readcount (readcount_w_children)
        "c": float,  # cumulative abundance (readcount_w_children)
        "a": float,  # species abundance
        "p": str,  # parent tax_id (absent for root)
        "%": float,  # percent mapped (sum_reads / n_mapped)
        "x": int,  # raw readcount
        "y": int,  # cumulative raw readcount
        "k": int,  # kmers (include_kmers=True only)
        "m": int,  # uniquemers (include_kmers=True only)
    },
    total=False,
)


class ClassificationSummary(TypedDict):
    n_reads: int
    n_mapped: int
    pct_mapped: float
    data: dict[str, TaxonEntry]
    stats: dict[str, int | float]
    filtered: bool
    has_abundances: bool
    hosts: list[str]


class ClassificationTableRow(TypedDict):
    tax_id: str
    name: str
    rank: str
    parent_tax_id: str | None
    readcount: int
    readcount_w_children: int
    unfiltered_readcount: int
    unfiltered_readcount_w_children: int
    abundance: float | None
    abundance_w_children: float | None


class ClassificationTable(TypedDict):
    n_reads: int
    host_tax_ids: list[str]
    table: list[ClassificationTableRow]


class PyTaxonomy:
    """Pure-Python taxonomy tree.

    Replaces the compiled ``taxonomy`` Rust extension for environments
    (e.g. pyodide) where native extensions are unavailable.
    """

    def __init__(self, nodes: dict[str, TaxonomyNode], children_map: dict[str, list[str]]) -> None:
        self._nodes = nodes
        self._children_map = children_map

    @classmethod
    def from_json(cls, json_str: str) -> PyTaxonomy:
        """Build a PyTaxonomy from a raw results JSON string.

        Parameters
        ----------
        json_str : str
            Full contents of a raw ``results.json`` file.  The ``taxonomy``
            key must contain a nested ``{id, name, rank, children: [...]}``
            tree.
        """
        data = json.loads(json_str)
        return cls.from_nested_tree(data["taxonomy"])

    @classmethod
    def from_nested_tree(cls, root: dict) -> PyTaxonomy:
        """Build a PyTaxonomy from a nested ``{id, name, rank, children}`` tree.

        Parameters
        ----------
        root : dict
            Root node of the nested taxonomy tree as found in raw
            ``results.json`` files.
        """
        nodes: dict[str, TaxonomyNode] = {}
        children_map: dict[str, list[str]] = {}

        def _traverse(node: dict, parent_id: str | None) -> None:
            tax_id = str(node["id"])
            nodes[tax_id] = TaxonomyNode(
                tax_id=tax_id,
                name=node["name"],
                rank=node["rank"],
                parent_tax_id=parent_id,
            )
            children_map[tax_id] = [str(c["id"]) for c in node.get("children", [])]
            for child in node.get("children", []):
                _traverse(child, tax_id)

        _traverse(root, None)
        return cls(nodes, children_map)

    @classmethod
    def from_table(cls, table: list[dict]) -> PyTaxonomy:
        """Build a PyTaxonomy from the flat ``table`` list in a processed results dict.

        Parameters
        ----------
        table : list[dict]
            The ``table`` value from a processed classification results dict.
            Each entry must have at least ``tax_id``, ``name``, ``rank``, and
            ``parent_tax_id`` keys.
        """
        nodes: dict[str, TaxonomyNode] = {}
        children_map: dict[str, list[str]] = {}

        for row in table:
            tax_id = row["tax_id"]
            nodes[tax_id] = TaxonomyNode(
                tax_id=tax_id,
                name=row["name"],
                rank=row["rank"],
                parent_tax_id=row.get("parent_tax_id"),
            )
            children_map.setdefault(tax_id, [])

        for tax_id, node in nodes.items():
            parent = node.parent_tax_id
            if parent is not None and parent in nodes:
                children_map.setdefault(parent, []).append(tax_id)

        return cls(nodes, children_map)

    def children(self, tax_id: str) -> list[str]:
        """Return the direct child tax IDs of *tax_id*."""
        return self._children_map.get(tax_id, [])

    def lineage(self, tax_id: str) -> list[str]:
        """Return the lineage from *tax_id* up to (and including) the root.

        Returns
        -------
        list[str]
            Ordered ``[tax_id, parent, grandparent, …, root]``.
        """
        result: list[str] = []
        current: str | None = tax_id
        seen: set[str] = set()
        while current is not None and current not in seen:
            seen.add(current)
            result.append(current)
            node = self._nodes.get(current)
            if node is None:
                break
            current = node.parent_tax_id
        return result

    def node(self, tax_id: str) -> TaxonomyNode | None:
        """Return the node for *tax_id*, or ``None`` if not found."""
        return self._nodes.get(tax_id)

    def __getitem__(self, tax_id: str) -> TaxonomyNode:
        return self._nodes[tax_id]

    def __contains__(self, tax_id: object) -> bool:
        return tax_id in self._nodes


def dfs_postorder_taxonomy(tree: PyTaxonomy, root_id: str = "1") -> Generator[str, None, None]:
    """Yield tax IDs in DFS post-order (children before parents).

    Parameters
    ----------
    tree : PyTaxonomy
    root_id : str, optional
        Tax ID of the root node.  Defaults to ``"1"`` (NCBI root).

    Yields
    ------
    str
        Tax IDs in post-order.
    """
    visited: set[str] = set()
    visited.add(root_id)
    stack = [(root_id, iter(tree.children(root_id)))]
    while stack:
        parent, children = stack[-1]
        try:
            child = next(children)
            if child not in visited:
                visited.add(child)
                stack.append((child, iter(tree.children(child))))
        except StopIteration:
            stack.pop()
            if stack:
                yield parent
    yield root_id


def shannon_diversity(proportions: Iterable[float], base: float = math.e) -> float:
    """Compute the Shannon alpha-diversity index from a sequence of proportions.

    Parameters
    ----------
    proportions : Iterable[float]
        Non-negative proportions.  Zero values are ignored.
    base : float, optional
        Logarithm base.  Defaults to *e* (natural log).

    Returns
    -------
    float
        Shannon H.  Returns ``0.0`` for empty or all-zero input.
    """
    return -sum(n * math.log(n, base) for n in proportions if n > 0)


def simpson_diversity(proportions: Iterable[float]) -> float:
    """Compute the Simpson diversity index (1 − D) from a sequence of proportions.

    Parameters
    ----------
    proportions : Iterable[float]
        Non-negative proportions.  Zero values are ignored.

    Returns
    -------
    float
        Simpson 1−D (Gini–Simpson coefficient).  Returns ``0.0`` for empty input.
    """
    return 1 - sum(p**2 for p in proportions if p > 0)


def analysis_stats(analysis_data: dict[str, TaxonEntry]) -> dict[str, int | float]:
    """Compute diversity statistics from the compact results data dict.

    Parameters
    ----------
    analysis_data : dict[str, TaxonEntry]
        The ``data`` dict as produced by :func:`summarize_analysis_json`.

    Returns
    -------
    dict[str, int | float]
        Diversity statistics keyed by metric name.
    """
    stats: dict[str, int | float] = {}

    abundances = [a.get("a", 0) for a in analysis_data.values() if a.get("a", 0) > 0]
    stats["species_richness_abund"] = len(abundances)
    stats["shannon_abund"] = shannon_diversity(abundances)
    stats["simpson_abund"] = simpson_diversity(abundances)

    for rank in ["species", "genus", "family", "order", "class", "phylum"]:
        readcounts = [a["t"] for a in analysis_data.values() if a["r"] == rank and a["t"] > 0]
        total = sum(readcounts)
        norm = [i / total for i in readcounts]
        stats["shannon_" + rank] = shannon_diversity(norm)
        stats["simpson_" + rank] = simpson_diversity(norm)

    return stats


def summarize_analysis_json(
    unprocessed_results: dict,
    include_kmers: bool = False,
) -> ClassificationSummary:
    """Summarise a raw One Codex classification results JSON.

    Transforms the raw ``results.json`` produced by the classification
    pipeline into the compact format returned by :func:`format_classification_table`.

    Parameters
    ----------
    unprocessed_results : dict
        Parsed raw ``results.json``.  Must contain ``n_reads``, ``n_mapped``,
        ``pct_mapped``, a ``results`` list, and a nested ``taxonomy`` tree.
    include_kmers : bool, optional
        If ``True``, include ``"k"`` (nKmers) and ``"m"`` (nUniquemers) fields
        in each :class:`TaxonEntry`.

    Returns
    -------
    ClassificationSummary
        A dict with ``n_reads``, ``n_mapped``, ``pct_mapped``, ``data``,
        ``stats``, ``filtered``, ``has_abundances``, and ``hosts`` keys.
    """
    node_data = {str(n["tax_id"]): n for n in unprocessed_results["results"]}

    json_taxa = json.dumps(unprocessed_results["taxonomy"])
    # Some analyses have a "synthetic" rank that the taxonomy library doesn't accept.
    json_taxa = json_taxa.replace('"rank": "synthetic"', '"rank": "no rank"')
    tree = PyTaxonomy.from_nested_tree(json.loads(json_taxa))

    data: dict[str, TaxonEntry] = {}

    nnls_ran = unprocessed_results.get("shares_db_dist", False) and unprocessed_results.get(
        "nnls_succeeded", False
    )
    otherwise_filtered = unprocessed_results.get("was_filtered", False)
    filtered = nnls_ran or otherwise_filtered

    n_mapped = (
        unprocessed_results["n_mapped"]
        if filtered
        else unprocessed_results.get("n_mapped_raw", unprocessed_results["n_mapped"])
    )

    if n_mapped > 0:
        for tax_id in dfs_postorder_taxonomy(tree):
            node = node_data.get(tax_id, {})
            child_tax_ids = tree.children(tax_id)

            # Fall back to readcount for old-style JSONs that lack rawReadcount.
            n_raw_reads = node.get("rawReadcount", node.get("readcount", 0))

            if filtered:
                n_reads = node.get("readcount", 0)
            else:
                n_reads = node.get("rawReadcount", node.get("readcount", 0))

            sum_reads = n_reads
            sum_raw_reads = n_raw_reads
            cumulative_abund = node.get("speciesAbundance", 0.0)

            if child_tax_ids:
                sum_reads += sum(data.get(t, {}).get("t", 0) for t in child_tax_ids)
                sum_raw_reads += sum(data.get(t, {}).get("y", 0) for t in child_tax_ids)
                cumulative_abund += sum(data.get(t, {}).get("c", 0) for t in child_tax_ids)

            if sum_reads == 0 and cumulative_abund == 0.0:
                continue

            entry: TaxonEntry = {
                "c": cumulative_abund,
                "t": sum_reads,
                "%": sum_reads / n_mapped,
                "i": tree[tax_id].name,
                "r": tree[tax_id].rank,
                "x": n_raw_reads,
                "y": sum_raw_reads,
            }

            if include_kmers:
                entry["k"] = node.get("nKmers", 0)
                entry["m"] = node.get("nUniquemers", 0)

            if n_reads > 0:
                entry["n"] = n_reads

            abundance = node.get("speciesAbundance", None)
            if abundance is not None:
                entry["a"] = abundance

            parent_ids = tree.lineage(tax_id)[1:]
            if parent_ids:
                entry["p"] = parent_ids[0]

            data[tax_id] = entry

    if "host_tax_ids" in unprocessed_results:
        hosts = [str(tax_id) for tax_id in unprocessed_results["host_tax_ids"]]
    else:
        hosts = ["9606"]

    return {
        "n_reads": unprocessed_results["n_reads"],
        "n_mapped": n_mapped,
        "pct_mapped": unprocessed_results["pct_mapped"] / 100,
        "data": data,
        "stats": analysis_stats(data),
        "filtered": filtered,
        "has_abundances": nnls_ran,
        "hosts": hosts,
    }


def format_classification_table(summary: ClassificationSummary) -> ClassificationTable:
    """Convert a compact classification summary to the flat table format served by the API.

    Parameters
    ----------
    summary : ClassificationSummary
        The compact dict produced by :func:`summarize_analysis_json`.

    Returns
    -------
    ClassificationTable
        A dict with ``n_reads``, ``host_tax_ids``, and ``table`` keys matching
        the ``/api/v1/classifications/{id}/results`` response format.
    """
    has_abundances = summary.get("has_abundances", False)
    table: list[ClassificationTableRow] = []

    for tax_id, entry in summary.get("data", {}).items():
        readcount = entry.get("n", 0)
        readcount_w_children = entry["t"]
        abundance = entry.get("a") if has_abundances else None
        abundance_w_children = entry.get("c") if has_abundances else None

        # Omit rows that carry no information.
        if (
            readcount < 1
            and readcount_w_children < 1
            and not abundance
            and not abundance_w_children
        ):
            continue

        table.append(
            {
                "tax_id": tax_id,
                "name": entry["i"],
                "rank": entry["r"],
                "parent_tax_id": entry.get("p"),
                "readcount": readcount,
                "readcount_w_children": readcount_w_children,
                "unfiltered_readcount": entry["x"],
                "unfiltered_readcount_w_children": entry["y"],
                "abundance": abundance,
                "abundance_w_children": abundance_w_children,
            }
        )

    return {
        "n_reads": summary["n_reads"],
        "host_tax_ids": summary["hosts"],
        "table": table,
    }
