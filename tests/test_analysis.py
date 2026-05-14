"""Tests for onecodex.analysis."""

import json
import math
import os

import pytest

from onecodex.classification_results import (
    RESULTS_FORMAT_VERSION,
    PyTaxonomy,
    format_classification_table,
    shannon_diversity,
    simpson_diversity,
    summarize_analysis_json,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "data", "results_unprocessed")


def _load(filename: str) -> dict:
    with open(os.path.join(DATA_DIR, filename)) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Test cases
# ---------------------------------------------------------------------------
# Each tuple: (present, n_reads, total_reads, parent_tax_id, abundance,
#              n_kmers, n_unique_kmers, raw_readcount, raw_readcount_w_children)

OLD_STYLE = {
    "file": "SRR1649162.10k.Se.ST1647.fastq.gz.results.json",
    "results": {
        "335795": (True, 1, 8, "186765", None, 0, 0, 1, 8),
        "59201": (True, 1266, 9327, "28901", None, 0, 0, 1266, 9327),
        "2": (True, 14, 9655, "131567", None, 0, 0, 14, 9655),
        "10239": (True, 0, 124, "1", None, 0, 0, 0, 124),
    },
}

NNLS_SUCCESS = {
    "file": "ERR589544.fastq.gz.results.json",
    "results": {
        "2": (True, 8030, 632693, "131567", None, 686801, 0, 8030, 633530),
        "28116": (True, 10568, 33684, "816", 0.1462, 449869, 242499, 10568, 33684),
        "997886": (True, 22594, 22594, "28116", None, 807857, 437140, 22594, 22594),
        "1410631": (False, None, None, None, None, 32, 24, None, None),
    },
}

NNLS_FAILED = {
    "file": "ERR747127.fastq.gz.results.json",
    "results": {
        "2": (True, 93062, 362215, "131567", None, 17758886, 0, 93062, 362215),
        "1309": (True, 44, 54, "1301", None, 793, 69, 44, 54),
        "1080071": (True, 21, 21, "1301", None, 138, 15, 21, 21),
    },
}

DN_SHARE_DB_DIST = {
    "file": "ERR589544_nnls_failed.fastq.gz.results.json",
    "results": {
        "2": (True, 8030, 695324, "131567", None, 686801, 0, 8030, 695324),
        "28116": (True, 10568, 33684, "816", None, 449869, 242499, 10568, 33684),
        "997886": (True, 22594, 22594, "28116", None, 807857, 437140, 22594, 22594),
        "1410631": (True, 4, 4, "186928", None, 32, 24, 4, 4),
    },
}

ALL_CASES = [OLD_STYLE, NNLS_SUCCESS, NNLS_FAILED, DN_SHARE_DB_DIST]


# ---------------------------------------------------------------------------
# RESULTS_FORMAT_VERSION
# ---------------------------------------------------------------------------


def test_results_format_version():
    assert RESULTS_FORMAT_VERSION == 9


# ---------------------------------------------------------------------------
# PyTaxonomy — from_nested_tree / from_json
# ---------------------------------------------------------------------------


@pytest.fixture
def old_style_json():
    return _load(OLD_STYLE["file"])


def test_taxonomy_from_nested_tree(old_style_json):
    tree = PyTaxonomy.from_nested_tree(old_style_json["taxonomy"])
    root = tree.node("1")
    assert root is not None
    assert root.name == "root"
    assert root.parent is None


def test_taxonomy_from_json(old_style_json):
    tree = PyTaxonomy.from_json(json.dumps(old_style_json))
    assert "1" in tree


# ---------------------------------------------------------------------------
# PyTaxonomy — from_table (flat format)
# ---------------------------------------------------------------------------

FLAT_FIXTURE_PATH = os.path.join(
    os.path.dirname(__file__),
    "data",
    "api",
    "v1",
    "classifications",
    "bef0bc57dd7f4c43",
    "results",
    "index.json",
)


@pytest.fixture
def flat_results():
    with open(FLAT_FIXTURE_PATH) as f:
        return json.load(f)


def test_taxonomy_from_table(flat_results):
    tree = PyTaxonomy.from_table(flat_results["table"])
    root = tree.node("1")
    assert root is not None
    assert root.name == "root"
    assert root.rank == "no rank"
    assert root.parent is None


def test_taxonomy_getitem(flat_results):
    tree = PyTaxonomy.from_table(flat_results["table"])
    node = tree["131567"]
    assert node.name == "cellular organisms"


def test_taxonomy_getitem_missing(flat_results):
    tree = PyTaxonomy.from_table(flat_results["table"])
    with pytest.raises(KeyError):
        _ = tree["99999999"]


def test_taxonomy_node_missing(flat_results):
    tree = PyTaxonomy.from_table(flat_results["table"])
    assert tree.node("99999999") is None


def test_taxonomy_children(flat_results):
    tree = PyTaxonomy.from_table(flat_results["table"])
    assert [n.id for n in tree.children("1")] == ["131567"]
    assert tree.children("1720194") == []


def test_taxonomy_lineage(flat_results):
    tree = PyTaxonomy.from_table(flat_results["table"])
    lineage = tree.lineage("1720194")
    assert lineage[0].id == "1720194"
    assert lineage[-1].id == "1"
    assert any(n.id == "1485" for n in lineage)


def test_taxonomy_contains(flat_results):
    tree = PyTaxonomy.from_table(flat_results["table"])
    assert "1" in tree
    assert "99999999" not in tree


# ---------------------------------------------------------------------------
# Diversity helpers
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "proportions, expected",
    [
        ([0.25, 0.25, 0.25, 0.25], math.log(4)),  # uniform -> H = ln(4)
        ([1.0], 0.0),  # single taxon -> H = 0
        ([], 0.0),  # empty
        ([0.8, 0.2], pytest.approx(0.5004, abs=1e-3)),
        ([0.1] * 10, pytest.approx(2.302, abs=1e-3)),
    ],
)
def test_shannon_diversity(proportions, expected):
    assert shannon_diversity(proportions) == expected


@pytest.mark.parametrize(
    "proportions, expected",
    [
        ([0.25, 0.25, 0.25, 0.25], 0.75),  # uniform -> 1 - 4*(0.25^2)
        ([1.0], 0.0),  # single taxon -> D = 0
        ([], 1.0),  # empty: 1 - sum([]) = 1 - 0 = 1
        ([0.8, 0.2], pytest.approx(0.320, abs=1e-3)),
        ([0.1] * 10, pytest.approx(0.900, abs=1e-3)),
    ],
)
def test_simpson_diversity(proportions, expected):
    assert simpson_diversity(proportions) == expected


# ---------------------------------------------------------------------------
# summarize_analysis_json — golden value tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "test_case", ALL_CASES, ids=["OLD_STYLE", "NNLS_SUCCESS", "NNLS_FAILED", "DN_SHARE_DB_DIST"]
)
def test_summarize_analysis_json_golden(test_case):
    ajson = _load(test_case["file"])
    summary = summarize_analysis_json(ajson)

    assert isinstance(summary["data"], dict)
    for key in summary["data"]:
        assert isinstance(key, str)

    expected_results = test_case["results"]
    for tax_id, expected in expected_results.items():
        present, n_reads, total_reads, parent, abundance, _, _, raw_rc, raw_rc_wc = expected

        if not present:
            assert tax_id not in summary["data"], f"tax_id {tax_id} should be absent"
            continue

        entry = summary["data"][tax_id]
        assert entry.get("n", 0) == n_reads, f"{tax_id}: n_reads"
        assert entry["t"] == total_reads, f"{tax_id}: total_reads"
        assert entry["p"] == parent, f"{tax_id}: parent"
        assert entry["x"] == raw_rc, f"{tax_id}: raw_readcount"
        assert entry["y"] == raw_rc_wc, f"{tax_id}: raw_readcount_w_children"

        if abundance is not None:
            assert entry["a"] == pytest.approx(abundance, abs=1e-4), f"{tax_id}: abundance"
            assert entry["c"] >= entry["a"]
        else:
            assert entry.get("a") is None
            assert entry["c"] >= 0.0


def test_summarize_analysis_json_nnls_abundances_sum_to_one():
    """When NNLS succeeded, species abundances must sum to ~1."""
    ajson = _load(NNLS_SUCCESS["file"])
    summary = summarize_analysis_json(ajson)

    assert summary["has_abundances"] is True
    total = sum(v.get("a", 0.0) for v in summary["data"].values())
    assert total == pytest.approx(1.0, abs=1e-4)


def test_summarize_analysis_json_filtered_flags():
    """Check filtered / has_abundances flags for each case."""
    cases = [
        (OLD_STYLE, False, False),
        (NNLS_SUCCESS, True, True),
        (NNLS_FAILED, False, False),
        (DN_SHARE_DB_DIST, False, False),
    ]
    for test_case, expected_filtered, expected_has_abundances in cases:
        ajson = _load(test_case["file"])
        summary = summarize_analysis_json(ajson)
        assert summary["filtered"] is expected_filtered, test_case["file"]
        assert summary["has_abundances"] is expected_has_abundances, test_case["file"]


# ---------------------------------------------------------------------------
# format_classification_table
# ---------------------------------------------------------------------------

EXPECTED_TABLE_ROW_KEYS = {
    "tax_id",
    "name",
    "rank",
    "parent_tax_id",
    "readcount",
    "readcount_w_children",
    "unfiltered_readcount",
    "unfiltered_readcount_w_children",
    "abundance",
    "abundance_w_children",
}


@pytest.mark.parametrize(
    "test_case", ALL_CASES, ids=["OLD_STYLE", "NNLS_SUCCESS", "NNLS_FAILED", "DN_SHARE_DB_DIST"]
)
def test_format_classification_table_shape(test_case):
    ajson = _load(test_case["file"])
    summary = summarize_analysis_json(ajson)
    result = format_classification_table(summary)

    assert result["n_reads"] == summary["n_reads"]
    assert result["host_tax_ids"] == summary["hosts"]
    assert isinstance(result["table"], list)
    assert len(result["table"]) > 0

    for row in result["table"]:
        assert set(row.keys()) == EXPECTED_TABLE_ROW_KEYS
        assert isinstance(row["tax_id"], str)
        assert row["readcount_w_children"] >= row["readcount"]
        assert row["unfiltered_readcount_w_children"] >= row["unfiltered_readcount"]


@pytest.mark.parametrize(
    "test_case", ALL_CASES, ids=["OLD_STYLE", "NNLS_SUCCESS", "NNLS_FAILED", "DN_SHARE_DB_DIST"]
)
def test_format_classification_table_golden(test_case):
    ajson = _load(test_case["file"])
    summary = summarize_analysis_json(ajson)
    result = format_classification_table(summary)

    by_tax_id = {row["tax_id"]: row for row in result["table"]}
    has_abundances = summary["has_abundances"]

    for tax_id, expected in test_case["results"].items():
        present, n_reads, total_reads, parent, abundance, _, _, raw_rc, raw_rc_wc = expected

        if not present:
            assert tax_id not in by_tax_id, f"tax_id {tax_id} should be absent from table"
            continue

        row = by_tax_id[tax_id]
        assert row["readcount"] == n_reads, f"{tax_id}: readcount"
        assert row["readcount_w_children"] == total_reads, f"{tax_id}: readcount_w_children"
        assert row["parent_tax_id"] == parent, f"{tax_id}: parent_tax_id"
        assert row["unfiltered_readcount"] == raw_rc, f"{tax_id}: unfiltered_readcount"
        assert (
            row["unfiltered_readcount_w_children"] == raw_rc_wc
        ), f"{tax_id}: unfiltered_readcount_w_children"

        if has_abundances and abundance is not None:
            assert row["abundance"] == pytest.approx(abundance, abs=1e-4), f"{tax_id}: abundance"
        else:
            assert row["abundance"] is None


def test_format_classification_table_no_zero_rows():
    """Rows where all counts are zero and no abundance should be omitted."""
    ajson = _load(NNLS_SUCCESS["file"])
    summary = summarize_analysis_json(ajson)
    result = format_classification_table(summary)

    for row in result["table"]:
        has_reads = row["readcount"] > 0 or row["readcount_w_children"] > 0
        has_abund = row["abundance"] or row["abundance_w_children"]
        assert has_reads or has_abund, f"zero-count row without abundance: {row['tax_id']}"


# ---------------------------------------------------------------------------
# analysis_stats
# ---------------------------------------------------------------------------


def test_analysis_stats_keys():
    ajson = _load(NNLS_SUCCESS["file"])
    summary = summarize_analysis_json(ajson)
    stats = summary["stats"]

    expected_keys = {
        "species_richness_abund",
        "shannon_abund",
        "simpson_abund",
        "shannon_species",
        "simpson_species",
        "shannon_genus",
        "simpson_genus",
        "shannon_family",
        "simpson_family",
        "shannon_order",
        "simpson_order",
        "shannon_class",
        "simpson_class",
        "shannon_phylum",
        "simpson_phylum",
    }
    assert set(stats.keys()) == expected_keys
