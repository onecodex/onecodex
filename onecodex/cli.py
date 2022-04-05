from __future__ import print_function
import click
import copy
import dateutil
from functools import partial
import logging
import os
import time
import warnings

from onecodex.api import Api
from onecodex.auth import _login, _logout, _remove_creds, login_required
from onecodex.lib.upload import DEFAULT_THREADS
from onecodex.metadata_upload import validate_appendables
from onecodex.scripts import subset_reads
from onecodex.utils import (
    click_path_autocomplete_helper,
    cli_resource_fetcher,
    CliLogFormatter,
    download_file_helper,
    valid_api_key,
    OPTION_HELP,
    progressbar,
    pprint,
    pretty_errors,
    run_via_threadpool,
    telemetry,
)
from onecodex.input_helpers import auto_detect_pairs, concatenate_multilane_files
from onecodex.version import __version__


# set the context for getting -h also
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
log = logging.getLogger("onecodex")


def warning_msg(message, category, filename, lineno, file=None, line=None):
    log.warning(message)


warnings.showwarning = warning_msg


# options
@click.group(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--api-key", callback=valid_api_key, help=OPTION_HELP["api_key"], metavar="<str:api_key>"
)
@click.option("--no-pretty-print", "no_pprint", is_flag=True, help=OPTION_HELP["no_pprint"])
@click.option("--verbose", "-v", is_flag=True, help=OPTION_HELP["verbose"])
@click.option(
    "--telemetry/--no-telemetry", is_flag=True, default=True, help=OPTION_HELP["telemetry"]
)
@click.version_option(version=__version__)
@click.pass_context
@telemetry
def onecodex(ctx, api_key, no_pprint, verbose, telemetry):
    """One Codex v1 API command line interface."""
    # Setup log formatter. TODO: Evaluate click-log instead
    log_formatter = CliLogFormatter()
    log.setLevel(logging.INFO)
    log.handlers[0].setFormatter(log_formatter)

    # set up the context for sub commands
    click.Context.get_usage = click.Context.get_help
    ctx.obj = {}
    ctx.obj["API_KEY"] = api_key
    ctx.obj["NOPPRINT"] = no_pprint
    ctx.obj["TELEMETRY"] = telemetry

    if verbose:
        log.setLevel(logging.DEBUG)


@onecodex.group("scripts", help="Assorted utility scripts")
def scripts():
    pass


scripts.add_command(subset_reads.cli, "subset_reads")

# TODO: remove filter_reads which is deprecated in favor of subset_reads
filter_reads = copy.deepcopy(subset_reads.cli)
filter_reads.hidden = True
scripts.add_command(filter_reads, "filter_reads")


@onecodex.group("documents", help="Access files in the Document Portal")
def documents():
    pass


@click.command("list", help="List files available to you")
@click.option("--json", is_flag=True, default=False, help="Output JSON instead of prettified table")
@click.pass_context
@telemetry
@login_required
def documents_list(ctx, json):
    docs_list = cli_resource_fetcher(ctx, "documents", [], print_results=json)

    if json:
        return

    if not docs_list:
        click.echo("You haven't uploaded any files yet, and no files have been shared with you.")
    else:

        def _size_formatter(size):
            suffix = "B"
            if size > 1e9:
                suffix = "GB"
                size /= 1e9
            elif size >= 1e6:
                suffix = "MB"
                size /= 1e6
            elif size >= 1e3:
                suffix = "KB"
                size /= 1e3

            return "%g %s" % (round(size, 2), suffix)

        formatters = ["%-18s", "%-34s", "%-25s", "%-11s", "%-12s"]
        table = [
            ["ID", "Name", "Owner", "Size", "Created On"],
            ["-" * 16, "-" * 32, "-" * 23, "-" * 9, "-" * 10],
        ]

        docs_list = sorted(
            docs_list,
            reverse=True,
            key=lambda x: time.mktime(dateutil.parser.parse(x["created_at"]).timetuple()),
        )

        for doc in docs_list:
            fname = doc["filename"]
            owner = doc["uploader"].email
            table.append(
                [
                    doc["$uri"].split("/")[-1],
                    fname if len(fname) <= 32 else fname[:29] + "...",
                    owner if len(owner) <= 23 else owner[:20] + "...",
                    _size_formatter(doc["size"]) if doc["size"] else "N/A",
                    dateutil.parser.parse(doc["created_at"]).strftime("%Y-%m-%d"),
                ]
            )

        for row in table:
            formatted_row = []
            for formatter, content in zip(formatters, row):
                formatted_row.append(formatter % content)
            click.echo("".join(formatted_row))


@click.command("upload", help="Upload a file to the Document Portal")
@click.option(
    "--max-threads",
    default=DEFAULT_THREADS,
    help=OPTION_HELP["max_threads"],
    metavar="<int:threads>",
)
@click.argument(
    "files",
    nargs=-1,
    required=False,
    type=click.Path(exists=True),
    shell_complete=partial(click_path_autocomplete_helper, directory=False),
)
@click.pass_context
@pretty_errors
@telemetry
@login_required
def documents_upload(ctx, max_threads, files):
    """Upload a document file (of any type) to One Codex."""
    if len(files) == 0:
        click.echo(ctx.get_help())
        return

    files = list(files)

    bar = click.progressbar(length=sum([os.path.getsize(x) for x in files]), label="Uploading... ")
    run_via_threadpool(
        ctx.obj["API"].Documents.upload,
        files,
        {"progressbar": bar},
        max_threads=8 if max_threads > 8 else max_threads,
        graceful_exit=False,
    )


@click.command("download", help="Download a file that has been shared with you")
@click.argument("file_id", nargs=1, required=True)
@click.option(
    "--output-document",
    "-O",
    "file_path",
    help="Write document to PATH",
    required=False,
    type=click.Path(dir_okay=False, allow_dash=True),
    shell_complete=partial(click_path_autocomplete_helper, filename=False),
)
@click.pass_context
@pretty_errors
@telemetry
@login_required
def documents_download(ctx, file_id, file_path):
    doc_obj = ctx.obj["API"].Documents.get(file_id)

    if not doc_obj:
        log.error("Could not find document {} (404 status code)".format(file_id))
        ctx.exit(1)

    if file_path == "-" or file_path == "/dev/stdout":
        doc_obj.download(file_obj=open("/dev/stdout", "wb"), progressbar=False)
    else:
        path = doc_obj.download(path=file_path, progressbar=True)
        click.echo("{} saved to {}".format(file_id, path), err=True)


documents.add_command(documents_list, "list")
documents.add_command(documents_upload, "upload")
documents.add_command(documents_download, "download")


# resources
@onecodex.command("analyses")
@click.argument("analyses", nargs=-1, required=False)
@click.pass_context
@telemetry
@login_required
def analyses(ctx, analyses):
    """Retrieve performed analyses."""
    cli_resource_fetcher(ctx, "analyses", analyses)


@onecodex.command("classifications")
@click.option("--read-level", "readlevel", is_flag=True, help=OPTION_HELP["readlevel"])
@click.option(
    "--read-level-path",
    "readlevel_path",
    type=click.Path(),
    default="./",
    help=OPTION_HELP["readlevel_path"],
    shell_complete=partial(click_path_autocomplete_helper, filename=False),
)
@click.option("--results", "results", is_flag=True, help=OPTION_HELP["results"])
@click.pass_context
@click.argument("classifications", nargs=-1, required=False)
@telemetry
@login_required
def classifications(ctx, classifications, results, readlevel, readlevel_path):
    """Retrieve performed metagenomic classifications."""

    # basic operation -- just print
    if not readlevel and not results:
        cli_resource_fetcher(ctx, "classifications", classifications)

    # fetch the results
    elif not readlevel and results:
        if len(classifications) != 1:
            log.error("Can only request results data on one Classification at a time")
        else:
            classification = ctx.obj["API"].Classifications.get(classifications[0])
            if not classification:
                log.error(
                    "Could not find classification {} (404 status code)".format(classifications[0])
                )
                return
            results = classification.results(json=True)
            pprint(results, ctx.obj["NOPPRINT"])

    # fetch the readlevel
    elif readlevel is not None and not results:
        if len(classifications) != 1:
            log.error("Can only request read-level data on one Classification at a time")
        else:
            classification = ctx.obj["API"].Classifications.get(classifications[0])
            if not classification:
                log.error(
                    "Could not find classification {} (404 status code)".format(classifications[0])
                )
                return
            tsv_url = classification._readlevel()["url"]
            log.info("Downloading tsv data from: {}".format(tsv_url))
            download_file_helper(tsv_url, readlevel_path)

    # both given -- complain
    else:
        log.error("Can only request one of read-level data or results data at a time")


@onecodex.command("panels")
@click.pass_context
@click.argument("panels", nargs=-1, required=False)
@telemetry
@login_required
def panels(ctx, panels):
    """Retrieve performed in silico panel results."""
    cli_resource_fetcher(ctx, "panels", panels)


@onecodex.command("samples")
@click.pass_context
@click.argument("samples", nargs=-1, required=False)
@telemetry
@login_required
def samples(ctx, samples):
    """Retrieve uploaded samples."""
    cli_resource_fetcher(ctx, "samples", samples)


# utilities
@onecodex.group("download", help="Download data from One Codex.")
def download_group():
    pass


@download_group.command("samples")
@click.argument(
    "outdir",
    type=click.Path(file_okay=False, dir_okay=True, writable=True),
    nargs=1,
    required=True,
)
@click.option(
    "--project", help="Filter to samples in a given project. Can be project name or UUID."
)
@click.option(
    "-t", "--tags", multiple=True, help="Filter to samples that include *any* of these tag names."
)
@click.option(
    "--prompt/--no-prompt",
    is_flag=True,
    default=True,
    help="Prompt for confirmation before downloading a large number of samples. Setting --no-prompt "
    "will allow running without any user intervention, e.g. in a script.",
)
@click.pass_context
@pretty_errors
@telemetry
@login_required
def download_samples_command(ctx, outdir, project, tags, prompt):
    """Download FASTA/Q files from One Codex.

    Samples may optionally be filtered by project and/or tags. By default, all samples in your
    account will be downloaded.

    OUTDIR is the name of the output directory where downloaded files will be saved.

    """
    from onecodex.lib.download import download_samples

    download_samples(
        ctx.obj["API"],
        outdir,
        project_name_or_id=project,
        tag_names=tags,
        prompt=prompt,
        progressbar=True,
    )


@onecodex.command("upload")
@click.argument(
    "files",
    nargs=-1,
    required=False,
    type=click.Path(exists=True),
    shell_complete=partial(click_path_autocomplete_helper, directory=False),
)
@click.option("--max-threads", default=4, help=OPTION_HELP["max_threads"], metavar="<int:threads>")
@click.option(
    "--coerce-ascii",
    is_flag=True,
    default=False,
    help="automatically rename unicode filenames to ASCII",
)
@click.option(
    "--forward",
    type=click.Path(exists=True),
    help=OPTION_HELP["forward"],
    shell_complete=partial(click_path_autocomplete_helper, directory=False),
)
@click.option(
    "--reverse",
    type=click.Path(exists=True),
    help=OPTION_HELP["reverse"],
    shell_complete=partial(click_path_autocomplete_helper, directory=False),
)
@click.option("--prompt/--no-prompt", is_flag=True, help=OPTION_HELP["prompt"], default=True)
@click.option("--tag", "-t", "tags", multiple=True, help=OPTION_HELP["tag"])
@click.option("--metadata", "-md", multiple=True, help=OPTION_HELP["metadata"])
@click.option("--project", "-p", "project_id", help=OPTION_HELP["project"])
@click.option("--sample-id", help=OPTION_HELP["sample_id"])
@click.option("--external-sample-id", help=OPTION_HELP["external_sample_id"])
@click.pass_context
@pretty_errors
@telemetry
@login_required
def upload(
    ctx,
    files,
    max_threads,
    coerce_ascii,
    forward,
    reverse,
    prompt,
    tags,
    metadata,
    project_id,
    sample_id,
    external_sample_id,
):
    """Upload a FASTA or FASTQ (optionally gzip'd) to One Codex."""
    appendables = {}
    if tags:
        appendables["tags"] = []
        for tag in tags:
            appendables["tags"].append(tag)

    if metadata:
        appendables["metadata"] = {}
        for metadata_kv in metadata:
            split_metadata = metadata_kv.split("=", 1)
            if len(split_metadata) > 1:
                metadata_value = split_metadata[1]
                appendables["metadata"][split_metadata[0]] = metadata_value

    appendables = validate_appendables(appendables, ctx.obj["API"])

    if (forward or reverse) and not (forward and reverse):
        click.echo("You must specify both forward and reverse files", err=True)
        ctx.exit(1)

    if forward and reverse:
        if len(files) > 0:
            click.echo(
                "You may not pass a FILES argument when using the "
                " --forward and --reverse options.",
                err=True,
            )
            ctx.exit(1)
        files = [(forward, reverse)]
    elif len(files) == 0:
        click.echo(ctx.get_help())
        return
    else:
        files_set = set(files)
        if files_set.symmetric_difference(files):
            click.echo(
                "Duplicate filenames detected in command line--please specific each file only once",
                err=True,
            )
            ctx.exit(1)

        files = auto_detect_pairs(files, prompt)

    files = concatenate_multilane_files(files, prompt)

    total_size = sum(
        [
            (os.path.getsize(x[0]) + os.path.getsize(x[1]))
            if isinstance(x, tuple)
            else os.path.getsize(x)
            for x in files
        ]
    )

    upload_kwargs = {
        "metadata": appendables["valid_metadata"],
        "tags": appendables["valid_tags"],
        "project": project_id,
        "coerce_ascii": coerce_ascii,
        "progressbar": progressbar(length=total_size, label="Uploading..."),
        "sample_id": sample_id,
        "external_sample_id": external_sample_id,
    }

    if (sample_id or external_sample_id) and len(files) > 1:
        click.echo(
            "Please only specify a single file or pair of files to upload if using `sample_id` or `external_sample_id`",
            err=True,
        )
        ctx.exit(1)

    run_via_threadpool(
        ctx.obj["API"].Samples.upload,
        files,
        upload_kwargs,
        max_threads=8 if max_threads > 8 else max_threads,
        graceful_exit=False,
    )


@onecodex.command("login")
@click.pass_context
@telemetry
def login(ctx):
    """Add an API key (saved in ~/.onecodex)."""
    base_url = os.environ.get("ONE_CODEX_API_BASE", "https://app.onecodex.com")
    if not ctx.obj["API_KEY"]:
        _login(base_url)
    else:
        email = _login(base_url, api_key=ctx.obj["API_KEY"])
        ocx = Api(api_key=ctx.obj["API_KEY"], telemetry=ctx.obj["TELEMETRY"])

        # TODO: This should be protected or built in as a first class resource
        # with, e.g., connection error catching (it's not part of our formally documeted API at the moment)
        if ocx._client.Account.instances()["email"] != email:
            click.echo("Your login credentials do not match the provided email!", err=True)
            _remove_creds()
            ctx.exit(1)


@onecodex.command("logout")
@click.pass_context
@telemetry
def logout(ctx):
    """Delete an API key (saved in ~/.onecodex)."""
    _logout()
