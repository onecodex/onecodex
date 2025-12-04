import json
import os
import csv
import pytest
from onecodex.utils import use_tempdir
from onecodex.scripts.export_functional_metric import (
    LongFunctionalResultExporter,
    WideFunctionalResultExporter,
)


def load_filtered_functional_results_jsons(taxa_stratified: bool = False) -> list[dict]:
    regular = [
        "functional_pathways_abundance.gz",
        "functional_pathways_abundance2.gz",
        "functional_pathways_abundance3.gz",
    ]
    stratified = [
        "functional_pathways_abundance_stratified.gz",
        "functional_pathways_abundance_stratified2.gz",
        "functional_pathways_abundance_stratified3.gz",
    ]

    files = stratified if taxa_stratified else regular
    return [
        json.loads(open(os.path.join("tests/data/files", filename)).read()) for filename in files
    ]


@pytest.mark.parametrize(
    "fmt,taxa_stratified,expected_header,expected_lines",
    [
        (
            "long",
            True,
            [
                "Functional Run ID",
                "Sample Name",
                "Sample ID",
                "Function ID",
                "Function Name",
                "Taxon ID",
                "Taxon Name",
                "Abundance",
            ],
            [
                [
                    "fr_0",
                    "Sample 0",
                    "s_0",
                    "PWY-6163",
                    "chorismate biosynthesis from 3-dehydroquinate",
                    "301302",
                    "Roseburia faecis",
                    "45.788",
                ],
                [
                    "fr_0",
                    "Sample 0",
                    "s_0",
                    "GLUDEG-I-PWY",
                    "GABA shunt",
                    "",
                    "unclassified",
                    "39.494",
                ],
                [
                    "fr_1",
                    "Sample 1",
                    "s_1",
                    "PWY-6123",
                    "inosine-5'-phosphate biosynthesis I",
                    "853",
                    "Faecalibacterium prausnitzii",
                    "44.667",
                ],
                [
                    "fr_1",
                    "Sample 1",
                    "s_1",
                    "PWY-6700",
                    "queuosine biosynthesis I (de novo)",
                    "46503",
                    "Parabacteroides merdae",
                    "22.375",
                ],
                [
                    "fr_2",
                    "Sample 2",
                    "s_2",
                    "COA-PWY",
                    "coenzyme A biosynthesis I (prokaryotic)",
                    "573",
                    "Klebsiella pneumoniae",
                    "1.568",
                ],
                [
                    "fr_2",
                    "Sample 2",
                    "s_2",
                    "P221-PWY",
                    "octane oxidation",
                    "573",
                    "Klebsiella pneumoniae",
                    "10.223",
                ],
            ],
        ),
        (
            "long",
            False,
            [
                "Functional Run ID",
                "Sample Name",
                "Sample ID",
                "Function ID",
                "Function Name",
                "Abundance",
            ],
            [
                ["fr_0", "Sample 0", "s_0", "ARO-PWY", "chorismate biosynthesis I", "709.749"],
                ["fr_0", "Sample 0", "s_0", "GLUDEG-I-PWY", "GABA shunt", "94.321"],
                ["fr_1", "Sample 1", "s_1", "PWY-6531", "mannitol cycle", "24.430"],
                ["fr_1", "Sample 1", "s_1", "PWY0-1479", "tRNA processing", "366.605"],
                ["fr_2", "Sample 2", "s_2", "PWY-6531", "mannitol cycle", "42.346"],
                ["fr_2", "Sample 2", "s_2", "PWY66-399", "gluconeogenesis III", "475.673"],
            ],
        ),
        (
            "wide",
            True,
            [
                "Function ID",
                "Function Name",
                "Taxon ID",
                "Taxon Name",
                "Sample 0 (run id fr_0)",
                "Sample 1 (run id fr_1)",
                "Sample 2 (run id fr_2)",
            ],
            [
                [
                    "PYRIDNUCSYN-PWY",
                    "NAD de novo biosynthesis I (from aspartate)",
                    "",
                    "unclassified",
                    "179.418",
                    "565.489",
                    "792.698",
                ],
                [
                    "NONOXIPENT-PWY",
                    "pentose phosphate pathway (non-oxidative branch) I",
                    "",
                    "Butyrivibrio crossotus CAG 259",
                    "",
                    "36.469",
                    "107.171",
                ],
                [
                    "PWY-6122",
                    "5-aminoimidazole ribonucleotide biosynthesis II",
                    "",
                    "Butyrivibrio crossotus CAG 259",
                    "",
                    "53.510",
                    "137.443",
                ],
                [
                    "PWY-6609",
                    "adenine and adenosine salvage III",
                    "",
                    "Oscillibacter sp 57 20",
                    "",
                    "20.396",
                    "101.820",
                ],
                [
                    "THISYNARA-PWY",
                    "superpathway of thiamine diphosphate biosynthesis III (eukaryotes)",
                    "",
                    "Butyrivibrio crossotus CAG 259",
                    "",
                    "29.121",
                    "73.934",
                ],
                [
                    "VALSYN-PWY",
                    "L-valine biosynthesis",
                    "",
                    "Butyrivibrio crossotus CAG 259",
                    "",
                    "116.516",
                    "187.354",
                ],
            ],
        ),
        (
            "wide",
            False,
            [
                "Function ID",
                "Function Name",
                "Sample 0 (run id fr_0)",
                "Sample 1 (run id fr_1)",
                "Sample 2 (run id fr_2)",
            ],
            [
                ["PWY-6531", "mannitol cycle", "34.531", "24.430", "42.346"],
                ["PWY-7234", "inosine-5'-phosphate biosynthesis III", "59.919", "", "537.451"],
                [
                    "RIBOSYN2-PWY",
                    "flavin biosynthesis I (bacteria and plants)",
                    "857.775",
                    "939.560",
                    "2486.907",
                ],
                [
                    "SER-GLYSYN-PWY",
                    "superpathway of L-serine and glycine biosynthesis I",
                    "362.450",
                    "455.298",
                    "530.553",
                ],
                [
                    "THISYNARA-PWY",
                    "superpathway of thiamine diphosphate biosynthesis III (eukaryotes)",
                    "374.351",
                    "293.187",
                    "323.551",
                ],
                [
                    "PWY0-1297",
                    "superpathway of purine deoxyribonucleosides degradation",
                    "",
                    "185.706",
                    "539.587",
                ],
            ],
        ),
    ],
)
def test_functional_long_exporter(fmt, taxa_stratified, expected_header, expected_lines):
    results = load_filtered_functional_results_jsons(taxa_stratified=taxa_stratified)
    with use_tempdir() as tempdir:
        out_file = os.path.join(tempdir, "functional_long_results.csv")

        exporter = (
            LongFunctionalResultExporter(
                out_path=out_file, taxa_stratified=taxa_stratified, metric="abundance"
            )
            if fmt == "long"
            else WideFunctionalResultExporter(
                out_path=out_file, taxa_stratified=taxa_stratified, metric="abundance"
            )
        )
        for idx, result in enumerate(results):
            exporter.consume_results(f"s_{idx}", f"Sample {idx}", f"fr_{idx}", result)
        exporter.produce_output()

        with open(out_file, "r", newline="") as tmp_file:
            reader = csv.reader(tmp_file)
            header = next(reader)
            assert header == expected_header

            lines = [line for line in reader]
            assert all(line in lines for line in expected_lines)
