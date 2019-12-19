import click
import bz2
import csv
import gzip
import io
import os
import warnings

from onecodex.auth import login_required
from onecodex.exceptions import OneCodexException, ValidationError
from onecodex.utils import download_file_helper, get_download_dest, pretty_errors


def fastfastq(file_path):
    filename, ext = os.path.splitext(file_path)

    if ext in {".gz", ".gzip"}:
        open_func = gzip.open
    elif ext in {".bz", ".bz2", ".bzip", ".bzip2"}:
        open_func = bz2.open
    else:
        open_func = io.open

    with open_func(file_path, "rb") as fp:
        idx = 0
        buf = []

        for line in fp:
            idx += 1
            buf.append(line)

            if len(buf) == 4:
                buf = b"".join(buf)

                if not buf.startswith(b"@"):
                    raise ValidationError("FASTQ record line {} does not start with @".format(idx))

                yield buf

                buf = []


def validating_parser(file_path, **io_kwargs):
    import skbio

    for rec in skbio.io.read(file_path, **io_kwargs):
        buf = io.BytesIO()
        rec.write(buf, **io_kwargs)
        buf.seek(0)
        yield buf.read()


def get_filtered_filename(file_path):
    filename = os.path.basename(file_path)
    filename, ext = os.path.splitext(filename)

    if ext in {".gz", ".gzip", ".bz", ".bzip", ".bz2", ".bzip2"}:
        filename, ext = os.path.splitext(filename)

    return "{}.filtered{}".format(filename, ext), ext


def make_taxonomy_dict(classification, parent=False):
    """Parse a Classification object into a `dict` mapping a tax_id to its children or parent.

    Notes
    -----
    Restricted to taxonomic IDs that are represented in the classification results.
    """

    tax_id_map = {}

    if parent:
        for row in classification.results()["table"]:
            if row["parent_tax_id"] is not None and row["tax_id"] is not None:
                tax_id_map[row["tax_id"]] = row["parent_tax_id"]
    else:
        for row in classification.results()["table"]:
            if row["parent_tax_id"] is not None and row["tax_id"] is not None:
                try:
                    tax_id_map[row["parent_tax_id"]].add(row["tax_id"])
                except KeyError:
                    tax_id_map[row["parent_tax_id"]] = set([row["tax_id"]])

    return tax_id_map


def recurse_taxonomy_map(tax_id_map, tax_id, parent=False):
    """Traverse up or down through a taxonomic tree.

    Finds all children (or parents) of the given tax_id.

    Parameters
    ----------
    tax_id_map : `dict`
        The output of make_taxonomy_map().
    tax_id : `int`
        The taxonomic ID of the node to begin the traversal from.
    parent : `bool`
        Return parents of the taxonomic node rather than children.
    """

    if parent:
        # TODO: allow filtering on tax_id and its parents, too
        raise NotImplementedError
    else:

        def _child_recurse(tax_id, visited):
            try:
                children = [tax_id] + list(tax_id_map[tax_id])
            except KeyError:
                children = [tax_id]

            for child in children:
                if child not in visited:
                    visited.append(child)
                    children.extend(_child_recurse(child, visited))

            return children

        return list(set(_child_recurse(tax_id, [])))


def too_many_fastx_records():
    raise ValidationError(
        "FASTX file(s) provided have more records than the classification results"
    )


@click.command(
    "subset_reads",
    help="Subset a FASTX file based on the taxonomic results from a CLASSIFICATION_ID. "
    "By default, reads (or pairs or reads) matching the given taxonomic ID are "
    "written to the path provided, ignoring low confidence taxonomic assignments. ",
)
@click.argument("classification_id")
@click.argument("fastx", type=click.Path())
@click.option(
    "-t",
    "--tax-id",
    "tax_ids",
    required=True,
    multiple=True,
    help="Subset reads mapping to tax IDs. May be passed multiple times.",
)
@click.option("-r", "--reverse", type=click.Path(), help="The reverse (R2) read file, optionally.")
@click.option(
    "--validate/--do-not-validate",
    is_flag=True,
    default=False,
    help="By default, skip validation of FASTQ records (which is slow). Enable if you are uncertain you are working with a valid FASTQ file",
)
@click.option(
    "--with-children",
    default=False,
    is_flag=True,
    help="Match child taxa of those given with -t (e.g., all strains of E. coli)",
)
@click.option(
    "--subset-pairs-independently",
    default=False,
    is_flag=True,
    help="By default, if either read in a pair matches, both will be retained in the subset "
    "file. With this option, R1 and R2 files will be evaluated independently. Note "
    "that the subset output FASTQs are *not* guaranteed to have the same number of "
    "reads!",
)
@click.option(
    "--exclude-reads",
    default=False,
    is_flag=True,
    help="By default, matching reads are kept. Choose this option to instead output "
    "reads that do *not* match.",
)
@click.option(
    "--include-lowconf",
    default=False,
    is_flag=True,
    help="By default, reads with low confidence taxonomic assignments are ignored. "
    "Choose this option to include them.",
)
@click.option(
    "-o", "--out", default=".", type=click.Path(), help="Where to save the filtered outputs"
)
@click.pass_context
@pretty_errors
@login_required
def cli(
    ctx,
    classification_id,
    fastx,
    reverse,
    tax_ids,
    with_children,
    subset_pairs_independently,
    exclude_reads,
    include_lowconf,
    out,
    validate,
):
    if ctx.info_name == "filter_reads":
        warnings.warn(
            "filter_reads will be removed in a future version. Please use subset_reads instead!"
        )

    if not len(tax_ids):
        raise OneCodexException("You must supply at least one tax ID")

    # fetch classification result object from API
    classification = ctx.obj["API"].Classifications.get(classification_id)
    if classification is None:
        raise ValidationError("Classification {} not found.".format(classification_id))

    # if with children, expand tax_ids by referring to the taxonomic tree
    if with_children:
        tax_id_map = make_taxonomy_dict(classification)

        new_tax_ids = []

        for t_id in tax_ids:
            new_tax_ids.extend(recurse_taxonomy_map(tax_id_map, t_id))

        tax_ids = new_tax_ids

    tax_ids = set(tax_ids)

    # pull the classification result TSV
    tsv_url = classification._readlevel()["url"]
    readlevel_path = get_download_dest("./", tsv_url)
    if not os.path.exists(readlevel_path):
        download_file_helper(tsv_url, "./")
    else:
        click.echo("Using cached read-level results: {}".format(readlevel_path), err=True)

    # count the number of rows in the TSV file
    with gzip.open(readlevel_path, "rt") as tsv:
        try:
            tsv_row_count = 0
            for _ in tsv:
                tsv_row_count += 1
            tsv_row_count -= 1  # discount header line
        except EOFError:
            click.echo(
                "\nWe encountered an error while processing the read "
                "level results. Please delete {} and try again.".format(readlevel_path),
                err=True,
            )
            raise

    if reverse:
        if tsv_row_count % 2 != 0:
            raise ValidationError(
                "Classification results cannot have odd number of records if using --reverse/-r"
            )

        tsv_row_count = int(tsv_row_count / 2.0)

    # determine the name of the output file(s)
    filtered_filename, ext = get_filtered_filename(fastx)
    filtered_filename = os.path.join(out, filtered_filename)
    if reverse:
        rev_filtered_filename = get_filtered_filename(reverse)[0]
        rev_filtered_filename = os.path.join(out, rev_filtered_filename)

    if ext in {".fa", ".fna", ".fasta"}:
        io_kwargs = {"format": "fasta"}
    elif ext in {".fq", ".fastq"}:
        io_kwargs = {"format": "fastq", "variant": "illumina1.8"}
    else:
        raise OneCodexException(
            "{}: extension must be one of .fa, .fna, .fasta, .fq, .fastq".format(fastx)
        )

    # do the actual filtering
    save_msg = "Saving subsetted reads: {}".format(filtered_filename)
    if reverse:
        save_msg += " and {}".format(rev_filtered_filename)
    click.echo(save_msg, err=True)

    # see mainline/#3513. we must set idx=0 here for cases where the fastx file is empty
    idx = 0

    with click.progressbar(length=tsv_row_count) as bar, gzip.open(readlevel_path, "rt") as tsv:
        reader = csv.DictReader(tsv, delimiter="\t")

        if reverse:
            if not validate and io_kwargs["format"] == "fastq":
                fwd_iter = fastfastq(fastx)
                rev_iter = fastfastq(reverse)
            else:
                fwd_iter = validating_parser(fastx, **io_kwargs)
                rev_iter = validating_parser(reverse, **io_kwargs)

            with io.open(filtered_filename, "wb") as out_file, io.open(
                rev_filtered_filename, "wb"
            ) as rev_out_file:  # noqa

                for idx, (fwd, rev) in enumerate(zip(fwd_iter, rev_iter)):
                    if idx == tsv_row_count:
                        too_many_fastx_records()
                    if idx % 1000 == 0:
                        bar.update(1000)
                    row = next(reader)  # necessary to do it this way for py2 compat
                    row2 = next(reader)

                    if subset_pairs_independently:
                        if include_lowconf:
                            if exclude_reads:
                                if row["Tax ID"] not in tax_ids:
                                    out_file.write(fwd)
                                if row2["Tax ID"] not in tax_ids:
                                    rev_out_file.write(rev)
                            else:
                                if row["Tax ID"] in tax_ids:
                                    out_file.write(fwd)
                                if row2["Tax ID"] in tax_ids:
                                    rev_out_file.write(rev)
                        else:
                            if exclude_reads:
                                if (
                                    row.get("Passed Filter", "T") == "T"
                                    and row["Tax ID"] not in tax_ids
                                ):
                                    out_file.write(fwd)
                                if (
                                    row2.get("Passed Filter", "T") == "T"
                                    and row2["Tax ID"] not in tax_ids
                                ):
                                    rev_out_file.write(rev)
                            else:
                                if (
                                    row.get("Passed Filter", "T") == "T"
                                    and row["Tax ID"] in tax_ids
                                ):
                                    out_file.write(fwd)
                                if (
                                    row2.get("Passed Filter", "T") == "T"
                                    and row2["Tax ID"] in tax_ids
                                ):
                                    rev_out_file.write(rev)
                    else:
                        if include_lowconf:
                            if exclude_reads:
                                if row["Tax ID"] not in tax_ids or row2["Tax ID"] not in tax_ids:
                                    out_file.write(fwd)
                                    rev_out_file.write(rev)
                            else:
                                if row["Tax ID"] in tax_ids or row2["Tax ID"] in tax_ids:
                                    out_file.write(fwd)
                                    rev_out_file.write(rev)
                        else:
                            if exclude_reads:
                                if (
                                    row.get("Passed Filter", "T") == "T"
                                    and row["Tax ID"] not in tax_ids
                                ) or (
                                    row2.get("Passed Filter", "T") == "T"
                                    and row2["Tax ID"] not in tax_ids
                                ):
                                    out_file.write(fwd)
                                    rev_out_file.write(rev)
                            else:
                                if (
                                    row.get("Passed Filter", "T") == "T"
                                    and row["Tax ID"] in tax_ids
                                ) or (
                                    row2.get("Passed Filter", "T") == "T"
                                    and row2["Tax ID"] in tax_ids
                                ):
                                    out_file.write(fwd)
                                    rev_out_file.write(rev)
        else:
            if not validate and io_kwargs["format"] == "fastq":
                fwd_iter = fastfastq(fastx)
            else:
                fwd_iter = validating_parser(fastx, **io_kwargs)

            with io.open(filtered_filename, "wb") as out_file:
                for idx, (fwd, row) in enumerate(zip(fwd_iter, reader)):
                    if idx == tsv_row_count:
                        too_many_fastx_records()
                    if idx % 1000 == 0:
                        bar.update(1000)
                    if include_lowconf:
                        if exclude_reads:
                            if row["Tax ID"] not in tax_ids:
                                out_file.write(fwd)
                        else:
                            if row["Tax ID"] in tax_ids:
                                out_file.write(fwd)
                    else:
                        if exclude_reads:
                            if (
                                row.get("Passed Filter", "T") == "T"
                                and row["Tax ID"] not in tax_ids
                            ):
                                out_file.write(fwd)
                        else:
                            if row.get("Passed Filter", "T") == "T" and row["Tax ID"] in tax_ids:
                                out_file.write(fwd)

        if idx < tsv_row_count - 1:  # 0-based idx, 1-based tsv_row_count
            raise ValidationError(
                "FASTX file(s) provided have fewer records than the classification results"
            )

        bar.finish()
