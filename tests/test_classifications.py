import pytest

pytest.importorskip("pandas")  # noqa

from onecodex.models import Classifications


# mocks an old results tree where parent taxa can be missing
@pytest.fixture
def sample_tree_old():
    yield [
        # this goes first because we had a bug that only occurred if a
        # species-level result came before its parent in the results dict.
        {
            "abundance": 0.3,
            "name": "species 4",
            "parent_tax_id": "2",
            "rank": "species",
            "readcount": 30,
            "readcount_w_children": 30,
            "tax_id": "4",
        },
        {
            "abundance": None,
            "name": "root",
            "parent_tax_id": None,
            "rank": "no rank",
            "readcount": 10,
            "readcount_w_children": 100,
            "tax_id": "1",
        },
        {
            "abundance": None,
            "name": "genus 2",
            "parent_tax_id": "1",
            "rank": "genus",
            "readcount": 0,
            "readcount_w_children": 30,
            "tax_id": "2",
        },
        {
            "abundance": None,
            "name": "genus 3",
            "parent_tax_id": "1",
            "rank": "genus",
            "readcount": 0,
            "readcount_w_children": 60,
            "tax_id": "3",
        },
        {
            "abundance": 0.1,
            "name": "species 5",
            "parent_tax_id": "100",  # not in tree
            "rank": "species",
            "readcount": 0,
            "readcount_w_children": 0,
            "tax_id": "5",
        },
        {
            "abundance": 0.2,
            "name": "species 6",
            "parent_tax_id": "3",
            "rank": "species",
            "readcount": 20,
            "readcount_w_children": 20,
            "tax_id": "6",
        },
        {
            "abundance": 0.4,
            "name": "species 6",
            "parent_tax_id": "3",  # not in tree
            "rank": "species",
            "readcount": 40,
            "readcount_w_children": 40,
            "tax_id": "7",
        },
    ]


# mocks a new, complete results tree
@pytest.fixture
def sample_tree_new():
    yield [
        {
            "abundance": None,
            "name": "root",
            "parent_tax_id": None,
            "rank": "no rank",
            "readcount": 10,
            "readcount_w_children": 100,
            "tax_id": "1",
        },
        {
            "abundance": None,
            "name": "genus 2",
            "parent_tax_id": "1",
            "rank": "genus",
            "readcount": 0,
            "readcount_w_children": 30,
            "tax_id": "2",
        },
        {
            "abundance": None,
            "name": "genus 3",
            "parent_tax_id": "1",
            "rank": "genus",
            "readcount": 0,
            "readcount_w_children": 60,
            "tax_id": "3",
        },
        {
            "abundance": 0.3,
            "name": "species 4",
            "parent_tax_id": "2",
            "rank": "species",
            "readcount": 20,
            "readcount_w_children": 30,
            "tax_id": "4",
        },
        {
            "abundance": 0.1,
            "name": "species 5",
            "parent_tax_id": "2",
            "rank": "species",
            "readcount": 10,
            "readcount_w_children": 0,
            "tax_id": "5",
        },
        {
            "abundance": 0.2,
            "name": "species 6",
            "parent_tax_id": "3",
            "rank": "species",
            "readcount": 20,
            "readcount_w_children": 20,
            "tax_id": "6",
        },
        {
            "abundance": 0.4,
            "name": "species 6",
            "parent_tax_id": "3",  # not in tree
            "rank": "species",
            "readcount": 40,
            "readcount_w_children": 40,
            "tax_id": "7",
        },
    ]


def test_abundance_rollups(sample_tree_old, sample_tree_new):
    table = {t["tax_id"]: t for t in Classifications._append_abundance_rollups(sample_tree_new)}

    all_genera_abund = sum(
        x["abundance_w_children"] for x in table.values() if x["rank"] == "genus"
    )
    all_species_abund = sum(
        x["abundance_w_children"] for x in table.values() if x["rank"] == "species"
    )

    assert table["5"]["abundance"] == 0.1
    assert table["2"]["abundance_w_children"] == 0.4

    assert table["1"]["abundance_w_children"] == 1.0
    assert all_genera_abund == 1.0
    assert all_species_abund == 1.0

    with pytest.warns(UserWarning, match="no assigned reads"):
        table = {t["tax_id"]: t for t in Classifications._append_abundance_rollups(sample_tree_old)}

    all_genera_abund = sum(
        x["abundance_w_children"] for x in table.values() if x["rank"] == "genus"
    )
    all_species_abund = sum(
        x["abundance_w_children"] for x in table.values() if x["rank"] == "species"
    )

    assert table["5"]["abundance"] is None
    assert (
        round(table["2"]["abundance_w_children"], 4) == 0.3333
    )  # old abundance / (sum of abundances)

    assert table["1"]["abundance_w_children"] == 1.0
    assert all_genera_abund == 1.0
    assert all_species_abund == 1.0
