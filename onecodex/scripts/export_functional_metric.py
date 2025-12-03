import click
import csv
import re
import os
import shutil
import tempfile

from abc import ABC, abstractmethod
from datetime import datetime
from onecodex.exceptions import OneCodexException
from onecodex.auth import login_required
from onecodex.utils import pretty_errors
from onecodex.lib.helpers import RateLimiter, hash_to_hex
from onecodex.lib.download import get_project
from onecodex.lib.enums import FunctionalAnnotationsMetric, FunctionalAnnotations
from typing import TextIO

API_BASE = "/api/v1/functional_profiles"
SPECIES_RE = re.compile(r"s__(\w+)", re.IGNORECASE | re.UNICODE)

FUNCTIONAL_ANNOTATIONS = {fa.name.lower() for fa in FunctionalAnnotations}
PATHWAY_METRICS = {
    fm.name.lower()
    for fm in FunctionalAnnotationsMetric.metrics_for_annotation(FunctionalAnnotations.Pathways)
}
OTHER_METRICS = {
    fm.name.lower()
    for fm in FunctionalAnnotationsMetric.metrics_for_annotation(FunctionalAnnotations.Go)
}


def parse_functional_taxon_name(taxon_name: str) -> str:
    match = SPECIES_RE.search(taxon_name)
    if match:
        return " ".join(x for x in match.group(1).split("_") if x)
    return taxon_name


class BaseFunctionalResultExporter(ABC):
    @abstractmethod
    def consume_results(self, sample_id: str, sample_name: str, analysis_id: str, results: dict):
        pass

    @abstractmethod
    def produce_output(self):
        pass


class LongFunctionalResultExporter(BaseFunctionalResultExporter):
    def __init__(
        self,
        *,
        out_path: str | None = None,
        file_obj: TextIO | None = None,
        taxa_stratified: bool,
        metric: str,
    ):
        if out_path is None and file_obj is None:
            raise OneCodexException("Provide one of out_path or file_obj")
        elif out_path is not None and file_obj is not None:
            raise OneCodexException("Provide only one of out_path or file_obj")

        self.fd = open(out_path, "w", newline="") if out_path is not None else file_obj
        self.csv_writer = csv.writer(self.fd)
        self.taxa_stratified = taxa_stratified
        self.metric = FunctionalAnnotationsMetric(metric)

        headers = (
            [
                "Functional Run ID",
                "Sample Name",
                "Sample ID",
                "Function ID",
                "Function Name",
                "Taxon ID",
                "Taxon Name",
            ]
            if taxa_stratified
            else [
                "Functional Run ID",
                "Sample Name",
                "Sample ID",
                "Function ID",
                "Function Name",
            ]
        )
        headers.append(self.metric.name)
        self.csv_writer.writerow(headers)

    def consume_results(self, sample_id: str, sample_name: str, analysis_id: str, results: dict):
        for item in results.get("table", []):
            row = [analysis_id, sample_name, sample_id, item["id"], item["name"]]
            if self.taxa_stratified:
                row.append(str(item.get("taxon_id") or ""))
                row.append(parse_functional_taxon_name(item.get("taxon_name") or ""))

            value = item["value"]
            row.append(f"{value:.3f}" if value is not None else "")
            self.csv_writer.writerow(row)

    def produce_output(self):
        self.fd.close()


class WideFunctionalResultExporter(BaseFunctionalResultExporter):
    def __init__(
        self,
        *,
        out_path: str | None = None,
        file_obj: TextIO | None = None,
        taxa_stratified: bool,
        metric: str,
    ):
        if out_path is None and file_obj is None:
            raise OneCodexException("Provide one of out_path or file_obj")
        elif out_path is not None and file_obj is not None:
            raise OneCodexException("Provide only one of out_path or file_obj")

        self.fd = open(out_path, "w", newline="") if out_path is not None else file_obj
        self.csv_writer = csv.writer(self.fd)
        self.taxa_stratified = taxa_stratified
        self.metric = FunctionalAnnotationsMetric(metric)

        self.headers = (
            [
                "Function ID",
                "Function Name",
                "Taxon ID",
                "Taxon Name",
            ]
            if taxa_stratified
            else [
                "Function ID",
                "Function Name",
            ]
        )

        self.functional_run_positions = {}
        self.function_names = {}
        self.temp_files = {}
        self.temp_dir = tempfile.mkdtemp(prefix="func_exp_")

    def consume_results(self, sample_id: str, sample_name: str, analysis_id: str, results: dict):
        if analysis_id not in self.functional_run_positions:
            self.functional_run_positions[analysis_id] = len(self.headers)
            self.headers.append(f"{sample_name} (run id {analysis_id})")

        for item in results.get("table", []):
            function_id = item["id"]
            if function_id not in self.function_names:
                self.function_names[function_id] = item["name"]

            # Split data up either by function_id or by (function_id, taxon_name, taxon_id)
            if self.taxa_stratified:
                taxon_id = str(item.get("taxon_id") or "")
                taxon_name = parse_functional_taxon_name(item.get("taxon_name") or "")
                key = (function_id, taxon_name, taxon_id)
            else:
                key = function_id

            if key not in self.temp_files:
                tmp_path = os.path.join(self.temp_dir, f"{hash_to_hex(key)}.csv")
                self.temp_files[key] = tmp_path
            else:
                tmp_path = self.temp_files[key]

            value = item["value"]

            # Save data to a temp file in long format
            with open(tmp_path, "a", newline="") as tmp_file:
                tmp_writer = csv.writer(tmp_file)
                tmp_writer.writerow([analysis_id, f"{value:.3f}" if value is not None else ""])

    def produce_output(self):
        click.echo("Pivoting to wide format table", err=True)
        self.csv_writer.writerow(self.headers)
        # Convert data in long format in temp files to wide format
        try:
            for key, tmp_path in self.temp_files.items():
                if self.taxa_stratified:
                    function_id, taxon_name, taxon_id = key
                    row = [
                        function_id,
                        self.function_names.get(function_id, ""),
                        taxon_id,
                        taxon_name,
                    ] + [""] * len(self.functional_run_positions)
                else:
                    function_id = key
                    row = [function_id, self.function_names.get(function_id, "")] + [""] * len(
                        self.functional_run_positions
                    )
                with open(tmp_path, "r", newline="") as tmp_file:
                    reader = csv.reader(tmp_file)
                    # Read long format temp file into wide row
                    for analysis_id, value in reader:
                        pos = self.functional_run_positions.get(analysis_id)
                        if pos is not None:
                            row[pos] = value
                self.csv_writer.writerow(row)
        finally:
            shutil.rmtree(self.temp_dir, ignore_errors=True)
            self.fd.close()


@click.command(
    "export_functional_metric",
    help="Export a single CSV file containing functional metric data for requested samples. "
    "Both long and wide formats are supported.",
)
@click.option(
    "-o",
    "--out",
    required=False,
    help="A CSV filename for where to save the output",
)
@click.option(
    "-a",
    "--annotation",
    required=True,
    help=f"Choose a functional annotation that you'd like to pull from the functional results, one of: {FUNCTIONAL_ANNOTATIONS}.",
)
@click.option(
    "-m",
    "--metric",
    required=True,
    help=f"Choose a functional metric that you'd like to pull from the functional results. Choose from {PATHWAY_METRICS} for pathway annotations, or from {OTHER_METRICS} for other annotations.",
)
@click.option(
    "--not-taxa-stratified",
    required=False,
    default=False,
    is_flag=True,
    help="Skip taxa stratification of the results.",
)
@click.option(
    "-p",
    "--project",
    required=False,
    help="Provide the project from which the samples should be exported. Can be project name or id. Provide either this option or explicit sample ids.",
)
@click.option(
    "-s",
    "--sample-ids",
    required=False,
    help="Provide a comma separated list of sample ids for the samples that should be exported. Provide either this option or a project.",
)
@click.option(
    "--table-format",
    required=False,
    default="long",
    help="Choose either the wide format or long (default).",
)
@click.pass_context
@pretty_errors
@login_required
def cli(ctx, out, annotation, metric, not_taxa_stratified, project, sample_ids, table_format):
    if not project and not sample_ids:
        raise OneCodexException("Either --project or --sample-ids needs to be provided")
    elif project and sample_ids:
        raise OneCodexException("Only one of --project or --sample-ids can be provided")

    metric = metric.lower()
    annotation = annotation.lower()
    if annotation not in FUNCTIONAL_ANNOTATIONS:
        raise OneCodexException(f"Invalid annotation: '{annotation}'")
    if metric not in FunctionalAnnotationsMetric.metrics_for_annotation(annotation):
        raise OneCodexException(
            f"Metric: '{metric}' is not a valid metric for annotation: '{annotation}'"
        )

    if not out:
        now = datetime.now()
        out = f"functional_{annotation}_{metric}_{now.strftime('%Y-%m-%d_%H-%M')}.csv"

    ocx = ctx.obj["API"]
    if project:
        click.echo(f"Fetching samples for project: {project}", err=True)
        ocx_project = get_project(ocx, project)
        samples = ocx.Samples.where(project=ocx_project)
    else:
        click.echo(f"Fetching samples: {sample_ids}", err=True)
        samples = ocx.Samples.where({"uuid": {"$in": sample_ids.split(",")}})

    if not samples:
        raise OneCodexException("No matching samples found")

    sample_ids = [s.id for s in samples]
    functional_runs = ocx.FunctionalProfiles.where(sample=sample_ids)

    limiter = RateLimiter(max_actions=10, period=1.0)

    exporter = (
        LongFunctionalResultExporter(
            out_path=out, taxa_stratified=not not_taxa_stratified, metric=metric
        )
        if table_format == "long"
        else WideFunctionalResultExporter(
            out_path=out, taxa_stratified=not not_taxa_stratified, metric=metric
        )
    )

    with click.progressbar(functional_runs, label="Downloading functional results") as bar:
        for fr in bar:
            limiter.acquire()

            data = fr._filtered_results(
                annotation=annotation, metric=metric, taxa_stratified=not not_taxa_stratified
            )
            sample_name = fr.sample.metadata.name or fr.sample.filename
            exporter.consume_results(fr.sample.id, sample_name, fr.id, data)

    exporter.produce_output()
    click.echo(f"Exported {len(functional_runs)} functional results", err=True)
