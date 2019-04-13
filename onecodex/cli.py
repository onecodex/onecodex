from __future__ import print_function
import click
import copy
import dateutil
import logging
import os
import re
import sys
import time
import warnings


from onecodex.api import Api
from onecodex.auth import _login, _logout, _remove_creds, login_required
from onecodex.lib.upload import DEFAULT_THREADS, _file_size
from onecodex.metadata_upload import validate_appendables
from onecodex.scripts import subset_reads
from onecodex.utils import (
    cli_resource_fetcher,
    CliLogFormatter,
    download_file_helper,
    valid_api_key,
    OPTION_HELP,
    progressbar,
    pprint,
    pretty_errors,
    run_via_threadpool,
    warn_if_insecure_platform,
    telemetry,
)
from onecodex.version import __version__


# set the context for getting -h also
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

# Modify the root logger directly
log = logging.getLogger()
log.setLevel(logging.INFO)

# Setup log formatter. TODO: Evaluate click-log instead
log_formatter = CliLogFormatter()
stream_handler = logging.StreamHandler(sys.stderr)
stream_handler.setFormatter(log_formatter)
log.addHandler(stream_handler)


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
    """One Codex v1 API command line interface"""
    # set up the context for sub commands
    click.Context.get_usage = click.Context.get_help
    ctx.obj = {}
    ctx.obj["API_KEY"] = api_key
    ctx.obj["NOPPRINT"] = no_pprint
    ctx.obj["TELEMETRY"] = telemetry

    if verbose:
        log.setLevel(logging.DEBUG)

    # handle checking insecure platform, we let upload command do it by itself
    if ctx.invoked_subcommand != "upload":
        warn_if_insecure_platform()


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
@click.argument("files", nargs=-1, required=False, type=click.Path(exists=True))
@click.pass_context
@pretty_errors
@telemetry
@login_required
def documents_upload(ctx, max_threads, files):
    """Upload a document file (of any type) to One Codex"""
    if len(files) == 0:
        click.echo(ctx.get_help())
        return

    files = list(files)

    bar = click.progressbar(length=sum([_file_size(x) for x in files]), label="Uploading... ")
    run_via_threadpool(
        ctx.obj["API"].Documents.upload,
        files,
        {"progressbar": bar},
        max_threads=max_threads,
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
    """Retrieve performed analyses"""
    cli_resource_fetcher(ctx, "analyses", analyses)


@onecodex.command("classifications")
@click.option("--read-level", "readlevel", is_flag=True, help=OPTION_HELP["readlevel"])
@click.option(
    "--read-level-path",
    "readlevel_path",
    type=click.Path(),
    default="./",
    help=OPTION_HELP["readlevel_path"],
)
@click.option("--results", "results", is_flag=True, help=OPTION_HELP["results"])
@click.pass_context
@click.argument("classifications", nargs=-1, required=False)
@telemetry
@login_required
def classifications(ctx, classifications, results, readlevel, readlevel_path):
    """Retrieve performed metagenomic classifications"""

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
    """Retrieve performed in silico panel results"""
    cli_resource_fetcher(ctx, "panels", panels)


@onecodex.command("samples")
@click.pass_context
@click.argument("samples", nargs=-1, required=False)
@telemetry
@login_required
def samples(ctx, samples):
    """Retrieve uploaded samples"""
    cli_resource_fetcher(ctx, "samples", samples)


# utilites
@onecodex.command("upload")
@click.option("--max-threads", default=4, help=OPTION_HELP["max_threads"], metavar="<int:threads>")
@click.argument("files", nargs=-1, required=False, type=click.Path(exists=True))
@click.option(
    "--coerce-ascii",
    is_flag=True,
    default=False,
    help="automatically rename unicode filenames to ASCII",
)
@click.option("--forward", type=click.Path(exists=True), help=OPTION_HELP["forward"])
@click.option("--reverse", type=click.Path(exists=True), help=OPTION_HELP["reverse"])
@click.option("--prompt/--no-prompt", is_flag=True, help=OPTION_HELP["prompt"], default=True)
@click.option("--tag", "-t", "tags", multiple=True, help=OPTION_HELP["tag"])
@click.option("--metadata", "-md", multiple=True, help=OPTION_HELP["metadata"])
@click.option("--project", "-p", "project_id", help=OPTION_HELP["project"])
@click.pass_context
@pretty_errors
@telemetry
@login_required
def upload(
    ctx, files, max_threads, prompt, forward, reverse, tags, metadata, project_id, coerce_ascii
):
    """Upload a FASTA or FASTQ (optionally gzip'd) to One Codex"""

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
        files = list(files)

        # "intelligently" find paired files and tuple them
        paired_files = []
        single_files = set(files)

        for filename in files:
            # convert "read 1" filenames into "read 2" and check that they exist; if they do
            # upload the files as a pair, autointerleaving them
            pair = re.sub("[._][Rr]1[._]", lambda x: x.group().replace("1", "2"), filename)

            # we don't necessary need the R2 to have been passed in; we infer it anyways
            if pair != filename and os.path.exists(pair):
                if not prompt and pair not in single_files:
                    # if we're not prompting, don't automatically pull in files
                    # not in the list the user passed in
                    continue

                paired_files.append((filename, pair))

                if pair in single_files:
                    single_files.remove(pair)

                single_files.remove(filename)

        auto_pair = True

        if prompt and len(paired_files) > 0:
            pair_list = ""
            for p in paired_files:
                pair_list += "\n  {}  &  {}".format(os.path.basename(p[0]), os.path.basename(p[1]))

            answer = click.confirm(
                "It appears there are paired files:{}\nInterleave them after upload?".format(
                    pair_list
                ),
                default="Y",
            )

            if not answer:
                auto_pair = False

        if auto_pair:
            files = paired_files + list(single_files)

    total_size = sum(
        [
            (_file_size(x[0], uncompressed=True) + _file_size(x[1], uncompressed=True))
            if isinstance(x, tuple)
            else _file_size(x, uncompressed=False)
            for x in files
        ]
    )

    upload_kwargs = {
        "metadata": appendables["valid_metadata"],
        "tags": appendables["valid_tags"],
        "project": project_id,
        "coerce_ascii": coerce_ascii,
        "progressbar": progressbar(length=total_size, label="Uploading..."),
    }

    run_via_threadpool(
        ctx.obj["API"].Samples.upload,
        files,
        upload_kwargs,
        max_threads=max_threads,
        graceful_exit=False,
    )


@onecodex.command("login")
@click.pass_context
@telemetry
def login(ctx):
    """Add an API key (saved in ~/.onecodex)"""
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
    """Delete your API key (saved in ~/.onecodex)"""
    _logout()
