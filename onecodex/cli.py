from __future__ import print_function

import copy
import json
import logging
import os
import time
import warnings
from functools import partial

import click

from onecodex.api import Api
from onecodex.auth import (
    _login,
    _logout,
    _remove_creds,
    login_required,
)
from onecodex.input_helpers import (
    auto_detect_pairs,
    concatenate_multilane_files,
    concatenate_ont_groups,
)
from onecodex.lib.upload import DEFAULT_THREADS
from onecodex.metadata_upload import validate_appendables
from onecodex.scripts import export_functional_metric, interleave, subset_reads
from onecodex.utils import (
    OPTION_HELP,
    CliLogFormatter,
    cli_resource_fetcher,
    click_path_autocomplete_helper,
    download_file_helper,
    pprint,
    pretty_errors,
    progressbar,
    run_via_threadpool,
    telemetry,
    use_tempdir,
    valid_api_key,
)
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
scripts.add_command(interleave.cli, "interleave")
scripts.add_command(export_functional_metric.cli, "export_functional_metric")

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
            key=lambda x: time.mktime(x.created_at.timetuple()),
        )
        # breakpoint()

        for doc in docs_list:
            fname = doc.filename
            owner = doc.uploader.email
            table.append(
                [
                    doc.field_uri.split("/")[-1],
                    fname if len(fname) <= 32 else fname[:29] + "...",
                    owner if len(owner) <= 23 else owner[:20] + "...",
                    _size_formatter(doc.size) if doc.size else "N/A",
                    doc.created_at.strftime("%Y-%m-%d"),
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


# assets


@onecodex.group("assets", help="Manage Assets.")
def assets():
    pass


@onecodex.command("upload")
@click.option(
    "--max-threads",
    default=DEFAULT_THREADS,
    help=OPTION_HELP["max_threads"],
    metavar="<int:threads>",
)
@click.option(
    "--name",
    nargs=1,
    required=False,
    help=OPTION_HELP["name"],
)
@click.argument(
    "file",
    nargs=1,
    required=True,
    type=click.Path(exists=True),
    shell_complete=partial(click_path_autocomplete_helper, directory=False),
)
@click.pass_context
@pretty_errors
@telemetry
@login_required
def assets_upload(ctx, max_threads, file, name):
    """Upload an asset to One Codex."""
    if len(file) == 0 or name == "":
        click.echo(ctx.get_help())
        return

    bar = click.progressbar(length=os.path.getsize(file), label="Uploading... ")
    run_via_threadpool(
        ctx.obj["API"].Assets.upload,
        [file],
        {"progressbar": bar, "name": name},
        max_threads=8 if max_threads > 8 else max_threads,
        graceful_exit=False,
    )


@click.command("list", help="List assets available to you")
@click.option("--json", is_flag=True, default=False, help="Output JSON instead of prettified table")
@click.pass_context
@telemetry
@login_required
def assets_list(ctx, json):
    assets_data = cli_resource_fetcher(ctx, "assets", [], print_results=json)
    if json:
        return

    if not assets_data:
        click.echo("You haven't uploaded any assets yet.")
        return

    formatters = ["%-18s", "%-34s", "%-12s", "%-12s"]
    table = [
        ["ID", "Filename", "Status", "Created On"],
        ["-" * 16, "-" * 32, "-" * 10, "-" * 10],
    ]

    assets_data = sorted(
        assets_data,
        reverse=True,
        key=lambda x: x.created_at.timestamp(),
    )

    for asset in assets_data:
        fname = asset.filename
        table.append(
            [
                asset.id,
                fname if len(fname) <= 32 else fname[:29] + "...",
                asset.status,
                asset.created_at.strftime("%Y-%m-%d"),
            ]
        )

    for row in table:
        formatted_row = []
        for formatter, content in zip(formatters, row):
            formatted_row.append(formatter % content)
        click.echo("".join(formatted_row))


assets.add_command(assets_upload, "upload")
assets.add_command(assets_list, "list")


# resources
class _AnalysesGroup(click.Group):
    """Group that preserves ``analyses [ID]...`` while still dispatching subcommands.

    Click's variadic ``nargs=-1`` argument on a group would otherwise eat the subcommand
    name. We split args manually: if any arg matches a registered subcommand, dispatch to
    it; otherwise treat the full arg list as analysis IDs for the group callback.
    """

    def parse_args(self, ctx, args):
        for i, a in enumerate(args):
            if a in self.commands:
                # Temporarily strip the variadic Argument so MultiCommand's normal
                # subcommand dispatch works; supply the captured IDs ourselves.
                ids = tuple(args[:i])
                variadic = [p for p in self.params if getattr(p, "name", None) == "analyses"]
                self.params = [p for p in self.params if p not in variadic]
                try:
                    rest = super().parse_args(ctx, args[i:])
                finally:
                    self.params.extend(variadic)
                ctx.params["analyses"] = ids
                return rest
        return click.Command.parse_args(self, ctx, args)


@onecodex.group("analyses", cls=_AnalysesGroup, invoke_without_command=True)
@click.argument("analyses", nargs=-1, required=False)
@click.pass_context
@telemetry
@login_required
def analyses(ctx, analyses):
    """Retrieve performed analyses."""
    if ctx.invoked_subcommand is None:
        cli_resource_fetcher(ctx, "analyses", analyses)


@click.command("await")
@click.argument("analysis_id", nargs=1, required=True)
@click.option(
    "--timeout",
    "timeout_seconds",
    type=float,
    default=None,
    help="Maximum number of seconds to wait. Waits indefinitely if not set.",
)
@click.option(
    "--initial-interval",
    "initial_interval_seconds",
    type=click.IntRange(min=5),
    default=5,
    show_default=True,
    help="Seconds between the first polls (minimum 5).",
)
@click.option(
    "--max-interval",
    "max_interval_seconds",
    type=click.IntRange(min=5),
    default=120,
    show_default=True,
    help="Upper bound (in seconds) on the polling interval after backoff (minimum 5).",
)
@click.pass_context
@pretty_errors
@telemetry
@login_required
def analyses_await(
    ctx, analysis_id, timeout_seconds, initial_interval_seconds, max_interval_seconds
):
    """Poll an analysis until it reaches a terminal state."""
    analysis = ctx.obj["API"].Analyses.get(analysis_id)
    if not analysis:
        raise click.ClickException(f"Could not find analysis {analysis_id} (404 status code)")

    try:
        analysis.await_completion(
            timeout=timeout_seconds,
            initial_interval=initial_interval_seconds,
            max_interval=max_interval_seconds,
        )
    except TimeoutError as exc:
        raise click.ClickException(str(exc))

    click.echo(f"Analysis {analysis.id} complete (success={analysis.success}).")
    if analysis.error_msg:
        click.echo(f"Error: {analysis.error_msg}", err=True)
    if analysis.success is False:
        ctx.exit(1)


onecodex.add_command(analyses_await, "await")
analyses.add_command(analyses_await, "await")


@click.command("logs")
@click.argument("analysis_id", nargs=1, required=True)
@click.option(
    "--tail",
    "tail",
    type=click.IntRange(min=1),
    default=1000,
    show_default=True,
    help="Fetch only the last N log lines.",
)
@click.pass_context
@pretty_errors
@telemetry
@login_required
def analyses_logs(ctx, analysis_id, tail):
    """Fetch the job run logs for an analysis."""
    analysis = ctx.obj["API"].Analyses.get(analysis_id)
    if not analysis:
        raise click.ClickException(f"Could not find analysis {analysis_id} (404 status code)")
    click.echo(analysis.logs(tail=tail), nl=False)


analyses.add_command(analyses_logs, "logs")


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


@onecodex.command("workflows")
@click.pass_context
@click.argument("workflows", nargs=-1, required=False)
@telemetry
@login_required
def workflows(ctx, workflows):
    """Retrieve performed workflow results."""
    cli_resource_fetcher(ctx, "workflows", workflows)


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

    with use_tempdir() as tempdir:
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

            # Detecting ONT groups comes first as otherwise part of ONT group could
            # be mistaken for a paired file
            files = concatenate_ont_groups(files, prompt, tempdir)
            files = auto_detect_pairs(files, prompt)

        files = concatenate_multilane_files(files, prompt, tempdir)

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
        resp = ocx._client.get(f"{ocx._base_url}/api/v1/account")
        resp.raise_for_status()
        if resp.json()["email"] != email:
            click.echo("Your login credentials do not match the provided email!", err=True)
            _remove_creds()
            ctx.exit(1)


@onecodex.command("logout")
@click.pass_context
@telemetry
def logout(ctx):
    """Delete an API key (saved in ~/.onecodex)."""
    _logout()


# jobs
@onecodex.group("jobs", help="One Codex platform Job management.")
def jobs_group():
    pass


@jobs_group.command("run")
@click.argument(
    "job_id",
    nargs=1,
    required=True,
)
@click.argument(
    "sample_id",
    nargs=1,
    required=True,
)
@click.option(
    "-a",
    "--arg",
    "args",
    multiple=True,
    help="Additional runtime arguments, example: -a min_quality=1 -a adapter=AGATC.",
)
@click.option(
    "-d",
    "--dependency-override",
    "dependency_overrides",
    multiple=True,
    help=(
        "Override a dependency with an existing analysis, example: "
        "-d <analysis_id> or -d <analysis_id>=<relative_download_path>."
    ),
)
@click.option(
    "--populate-default-arguments/--no-populate-default-arguments",
    default=True,
    help="Populate job arguments with their defaults on the server (default: enabled).",
)
@click.option(
    "--await",
    "await_completion",
    is_flag=True,
    default=False,
    help="Block until the new analysis reaches a terminal state.",
)
@click.pass_context
@pretty_errors
@telemetry
@login_required
def jobs_run(
    ctx,
    job_id,
    sample_id,
    args,
    dependency_overrides,
    populate_default_arguments,
    await_completion,
):
    """Run a OneCodex job with optional arguments."""
    from onecodex.models.misc import DependencyOverride

    parsed_args = {}
    for arg in args:
        if "=" not in arg:
            raise click.BadParameter(
                f"Expected key=value format, got {arg!r}.", param_hint="-a/--arg"
            )
        key, value = arg.split("=", 1)
        parsed_args[key] = value

    parsed_dependencies = []
    for dep in dependency_overrides:
        analysis_id, _, download_path = dep.partition("=")
        if not analysis_id:
            raise click.BadParameter(
                f"Expected <analysis_id> or <analysis_id>=<download_path>, got {dep!r}.",
                param_hint="-d/--dependency-override",
            )
        dep_analysis = ctx.obj["API"].Analyses.get(analysis_id)
        parsed_dependencies.append(
            DependencyOverride(analysis=dep_analysis, download_path=download_path or None)
        )

    job = ctx.obj["API"].Jobs.get(job_id)
    sample = ctx.obj["API"].Samples.get(sample_id)

    run = job.run(
        sample,
        parsed_args,
        dependency_overrides=parsed_dependencies or None,
        populate_default_arguments=populate_default_arguments,
    )
    click.echo(f"Job run created successfully. New analysis ID: {run.id}")

    if await_completion:
        run.await_completion()
        click.echo(f"Analysis {run.id} complete (success={run.success}).")
        if run.error_msg:
            click.echo(f"Error: {run.error_msg}", err=True)
        if run.success is False:
            ctx.exit(1)
    else:
        click.echo(f"Get its status using `onecodex analyses {run.id}`", err=True)


JOB_TYPE_CHOICES = ["shell_script", "nextflow"]


def _parse_dependency_specs(deps: tuple[str, ...]) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    for dep in deps:
        job_id, sep, output_dir = dep.partition("=")
        if not sep or not job_id or not output_dir:
            raise click.BadParameter(
                f"Expected <job_id>=<output_dir>, got {dep!r}.",
                param_hint="-d/--dependency",
            )
        pairs.append((job_id, output_dir))
    return pairs


def _build_repository(url: str | None, tag: str | None) -> dict[str, str | None] | None:
    if url is None and tag is None:
        return None
    if url is None:
        raise click.BadParameter(
            "--repository-tag requires --repository-url.",
            param_hint="--repository-tag",
        )
    return {"url": url, "tag": tag}


def _job_kwargs_from_options(
    api,
    name,
    script_path,
    image_uri,
    job_type,
    description,
    cpu,
    ram_gb,
    storage_gb,
    inject_bearer_token,
    repository_url,
    repository_tag,
    assets,
    dependencies,
    arguments_schema_path,
    autorun_on_org_sample_upload=None,
):
    asset_objs = [api.Assets.get(aid) for aid in assets] if assets else None
    dep_objs = (
        [
            {"job": api.Jobs.get(jid), "output_dir": out}
            for jid, out in _parse_dependency_specs(dependencies)
        ]
        if dependencies
        else None
    )
    script = open(script_path).read() if script_path else None
    arguments_schema = json.load(open(arguments_schema_path)) if arguments_schema_path else None
    payload = {
        "name": name,
        "script": script,
        "image_uri": image_uri,
        "job_type": job_type,
        "description": description,
        "cpu": cpu,
        "ram_gb": ram_gb,
        "storage_gb": storage_gb,
        "inject_bearer_token": inject_bearer_token,
        "repository": _build_repository(repository_url, repository_tag),
        "assets": asset_objs,
        "dependencies": dep_objs,
        "arguments_schema": arguments_schema,
        "autorun_on_org_sample_upload": autorun_on_org_sample_upload,
    }
    return {k: v for k, v in payload.items() if v is not None}


_JOB_OPTIONS_COMMON = [
    click.option("--image-uri", default=None, help="Fully qualified OCI image URI."),
    click.option(
        "--job-type",
        type=click.Choice(JOB_TYPE_CHOICES),
        default=None,
        help="Type of job (shell_script or nextflow).",
    ),
    click.option("--description", default=None, help="Human-readable description (markdown)."),
    click.option("--cpu", type=float, default=None, help="CPU cores (shell_script jobs)."),
    click.option("--ram-gb", type=float, default=None, help="RAM in GB (shell_script jobs)."),
    click.option(
        "--storage-gb", type=float, default=None, help="Storage in GB (shell_script jobs)."
    ),
    click.option(
        "--inject-bearer-token/--no-inject-bearer-token",
        "inject_bearer_token",
        default=None,
        help="Inject the user's bearer token into the job environment.",
    ),
    click.option("--repository-url", default=None, help="HTTPS URL of the git repository."),
    click.option("--repository-tag", default=None, help="Git tag for the repository."),
    click.option(
        "--asset-id",
        "assets",
        multiple=True,
        metavar="ASSET_ID",
        help="Asset to associate with the job (repeatable).",
    ),
    click.option(
        "--dependency",
        "-d",
        "dependencies",
        multiple=True,
        metavar="JOB_ID=OUTPUT_DIR",
        help="Job dependency, e.g. -d <job_id>=<output_dir>. Repeatable.",
    ),
    click.option(
        "--arguments-schema",
        "arguments_schema_path",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help="Path to a JSON file containing the arguments schema (list of FieldGroup).",
    ),
]


def _apply_options(options):
    """Apply a shared list of click.option() decorators to a command.

    Reverses the list so help-text order matches the list's declaration order
    (Click stacks decorators bottom-up).
    """

    def decorator(func):
        for option in reversed(options):
            func = option(func)
        return func

    return decorator


@jobs_group.command("create")
@click.option("--name", required=True, help="Human-readable name of the job.")
@click.option(
    "--script",
    "script_path",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Path to the script file to be executed.",
)
@_apply_options(_JOB_OPTIONS_COMMON)
@click.pass_context
@pretty_errors
@telemetry
@login_required
def jobs_create(
    ctx,
    name,
    script_path,
    image_uri,
    job_type,
    description,
    cpu,
    ram_gb,
    storage_gb,
    inject_bearer_token,
    repository_url,
    repository_tag,
    assets,
    dependencies,
    arguments_schema_path,
):
    """Create a new Job."""
    if image_uri is None:
        raise click.BadParameter("--image-uri is required.", param_hint="--image-uri")

    kwargs = _job_kwargs_from_options(
        api=ctx.obj["API"],
        name=name,
        script_path=script_path,
        image_uri=image_uri,
        job_type=job_type,
        description=description,
        cpu=cpu,
        ram_gb=ram_gb,
        storage_gb=storage_gb,
        inject_bearer_token=inject_bearer_token,
        repository_url=repository_url,
        repository_tag=repository_tag,
        assets=assets,
        dependencies=dependencies,
        arguments_schema_path=arguments_schema_path,
    )
    job = ctx.obj["API"].Jobs.create(**kwargs)
    click.echo(f"Created job {job.id}: {job.name}")


@jobs_group.command("update")
@click.argument("job_id", nargs=1, required=True)
@click.option("--name", default=None, help="Human-readable name of the job.")
@click.option(
    "--script",
    "script_path",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Path to the script file to be executed.",
)
@_apply_options(_JOB_OPTIONS_COMMON)
@click.option(
    "--autorun-on-org-sample-upload/--no-autorun-on-org-sample-upload",
    "autorun_on_org_sample_upload",
    default=None,
    help="Auto-run this job on every newly uploaded org sample (admins only).",
)
@click.pass_context
@pretty_errors
@telemetry
@login_required
def jobs_update(
    ctx,
    job_id,
    name,
    script_path,
    image_uri,
    job_type,
    description,
    cpu,
    ram_gb,
    storage_gb,
    inject_bearer_token,
    repository_url,
    repository_tag,
    assets,
    dependencies,
    arguments_schema_path,
    autorun_on_org_sample_upload,
):
    """Update an existing Job."""
    if job_type is not None:
        # The server-side UpdateJobSchema does not accept job_type.
        raise click.BadParameter(
            "--job-type cannot be changed after job creation.", param_hint="--job-type"
        )

    kwargs = _job_kwargs_from_options(
        api=ctx.obj["API"],
        name=name,
        script_path=script_path,
        image_uri=image_uri,
        job_type=None,
        description=description,
        cpu=cpu,
        ram_gb=ram_gb,
        storage_gb=storage_gb,
        inject_bearer_token=inject_bearer_token,
        repository_url=repository_url,
        repository_tag=repository_tag,
        assets=assets,
        dependencies=dependencies,
        arguments_schema_path=arguments_schema_path,
        autorun_on_org_sample_upload=autorun_on_org_sample_upload,
    )

    if not kwargs:
        raise click.UsageError("No fields provided to update.")

    job = ctx.obj["API"].Jobs.get(job_id)
    job.update(**kwargs)
    click.echo(f"Updated job {job.id}: {job.name}")
