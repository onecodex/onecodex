import base64
import click
import concurrent.futures
from functools import partial, wraps
import json
import logging
import os
import platform
import re
import requests
import sys
import sentry_sdk

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO  # noqa

try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse

from onecodex.vendored.potion_client.converter import PotionJSONEncoder
from onecodex.exceptions import OneCodexException, UploadException
from onecodex.version import __version__


log = logging.getLogger("onecodex")

OPTION_HELP = {
    "api_key": "Manually provide a One Codex API key",
    "no_pprint": "Do not pretty-print JSON responses",
    "threads": "Do not use multiple background threads to upload files",  # noqa
    "max_threads": "Specify a different max # of upload threads (defaults to 4)",  # noqa
    "verbose": "Log extra information to STDERR",
    "results": "Get a JSON array of the metagenomic classification results table",
    "readlevel": "Get the read-level data as a .tsv file",
    "readlevel_path": (
        "Output path or directory for a .tsv file with the raw "
        "read-level results. Defaults to original filename in "
        "the current working directory"
    ),
    "clean": (
        "Automatically clean up FASTX records during upload. This removes tabs from "
        "headers and converts U to T. If this isn't passed, these will cause errors."
    ),
    "interleave": (
        "Do not automatically interleave paired end files. Note this normally happens "
        "during the upload, does not require additional disk space, and does not change "
        "the files themselves."
    ),
    "prompt": (
        "Manually prompt about automatic paired file interleaving. Setting --no-prompt "
        "will allow running without any user intervention, e.g. in a script."
    ),
    "validate": (
        "Do not validate the FASTA/Q file before uploading. Incompatible with automatic "
        "paired end interleaving (NOT RECOMMENDED)."
    ),
    "telemetry": "Send errors to One Codex?",
    "forward": "Specify a forward reads file",
    "reverse": "Specify a reverse reads file",
    "tag": (
        "NOTE: We recommend invoking the upload command separately for each "
        "sample to add sample-specific tags.\n\nAdd one or more tags to "
        'all uploaded samples e.g., `onecodex upload --tag "Cohort A" $FILE`.'
    ),
    "metadata": (
        "NOTE: We recommend invoking the upload command separately for each "
        "sample to add sample-specific metadata.\n\nAdd one or more "
        "metadata attributes to all uploaded samples, "
        "e.g. `onecodex upload --metadata starred=true --metadata "
        'platform="Illumina MiSeq" $FILE`. '
    ),
    "project": "Provide the name, short name, or alphanumeric UUID of a "
    "project to automatically add the samples to that project on "
    "upload. NOTE: does not currently support adding a sample to "
    "a public project. Projects are searched by UUID, then name, "
    "then short name in that order.",
    "sample_id": "Provide an ID for a sample that was previously 'pre-uploaded' along with metadata.",
    "external_sample_id": "Provide an external sample ID for a sample that was previously 'pre-uploaded' along with metadata.",
}

SUPPORTED_EXTENSIONS = [
    "fa",
    "fasta",
    "fq",
    "fastq",
    "fa.gz",
    "fasta.gz",
    "fq.gz",
    "fastq.gz",
    "fa.gzip",
    "fasta.gzip",
    "fq.gzip",
    "fastq.gzip",
]


def valid_api_key(ctx, param, value):
    """Ensure an API has valid length (this is a click callback)."""
    if value is not None and len(value) != 32:
        raise click.BadParameter(
            "API Key must be 32 characters long, not {}".format(str(len(value)))
        )
    else:
        return value


def pprint(j, no_pretty):
    """Print as formatted JSON."""
    if not no_pretty:
        click.echo(
            json.dumps(j, cls=PotionJSONEncoder, sort_keys=True, indent=4, separators=(",", ": "))
        )
    else:
        click.echo(j)


class CliLogFormatter(logging.Formatter):
    formats = {logging.DEBUG: "DEBUG: %(module)s: %(lineno)d: %(msg)s", logging.INFO: "\n%(msg)s"}

    def __init__(self):
        # Note: only style="%" is supported in Python 2 and so we cannot pass style= here
        super(CliLogFormatter, self).__init__(fmt="\n%(levelname)s: %(msg)s", datefmt=None)
        self.original_format = self._style._fmt if hasattr(self, "_style") else self._fmt

    def format(self, record):
        if hasattr(self, "_style"):
            # Python 3
            self._style._fmt = self.formats.get(record.levelno, self.original_format)
        else:
            # Python 3
            self._fmt = self.formats.get(record.levelno, self.original_format)
        result = logging.Formatter.format(self, record)
        return result


def cli_resource_fetcher(ctx, resource, uris, print_results=True):
    try:
        # analyses is passed, want Analyses
        resource_name = resource[0].upper() + resource[1:]
        if len(uris) == 0:

            # if non given fetch all
            log.debug("No %s IDs given, fetching all...", resource_name)
            instances = getattr(ctx.obj["API"], resource_name).all()
            log.debug("Fetched %i %ss", len(instances), resource)
            objs_to_return = [x._resource._properties for x in instances]
        else:
            uris = list(set(uris))
            log.debug("Fetching %s: %s", resource_name, ",".join(uris))

            instances = []
            for uri in uris:
                try:
                    instance = getattr(ctx.obj["API"], resource_name).get(uri)
                    if instance is not None:
                        instances.append(instance._resource._properties)
                    else:
                        log.error(
                            "Could not find {} {} (404 status code)".format(resource_name, uri)
                        )
                except requests.exceptions.HTTPError as e:
                    log.error(
                        "Could not find {} {} ({} status code)".format(
                            resource_name, uri, e.response.status_code
                        )
                    )
            objs_to_return = instances

        if print_results:
            pprint(objs_to_return, ctx.obj["NOPPRINT"])
        else:
            return objs_to_return
    except requests.exceptions.HTTPError:
        click.echo(
            "Failed to authenticate. Please check your API key "
            "or trying logging out and back in with `onecodex logout` "
            "and `onecodex login`."
        )


def get_download_dest(input_path, url):
    original_filename = urlparse(url).path.split("/")[-1]
    if os.path.isdir(input_path):
        local_full_path = os.path.join(input_path, original_filename)
    else:
        local_full_path = input_path
    return local_full_path


def download_file_helper(url, input_path):
    """Manage the chunked downloading of a file given an url."""
    r = requests.get(url, stream=True)
    if r.status_code != 200:
        log.error("Failed to download file: %s" % r.json()["message"])
    local_full_path = get_download_dest(input_path, r.url)
    original_filename = os.path.split(local_full_path)[-1]
    with open(local_full_path, "wb") as f:
        click.echo("Downloading {}".format(original_filename), err=True)
        for chunk in r.iter_content(chunk_size=1024):
            if chunk:  # filter out keep-alive new chunks
                f.write(chunk)
                f.flush()
    pprint("Successfully downloaded %s to %s" % (original_filename, local_full_path), True)  # noqa


def check_for_allowed_file(f):
    """Check a file extension against a list of supporting sequence file extensions."""
    for ext in SUPPORTED_EXTENSIONS:
        if f.endswith(ext):
            return True
    log.error("Failed upload: Not an allowed file extension: %s", f)
    raise SystemExit


def collapse_user(fp):
    """Convert a path back to ~/ from expanduser()."""
    home_dir = os.path.expanduser("~")
    abs_path = os.path.abspath(fp)
    return abs_path.replace(home_dir, "~")


def _setup_sentry_for_ipython():
    from IPython import get_ipython

    ip = get_ipython()

    def custom_exc(shell, etype, evalue, tb, tb_offset=None):
        # Show original error
        shell.showtraceback((etype, evalue, tb), tb_offset=tb_offset)

        # Then send to Sentry
        sentry_sdk.capture_exception()

    if ip is not None:
        ip.set_custom_exc((Exception,), custom_exc)
        return True

    return False


def init_sentry(user_context=None, extra_context=None):
    if os.environ.get("ONE_CODEX_NO_TELEMETRY") is None:
        default_dsn = base64.b64decode(
            b"aHR0cHM6Ly9mYTU1ODliNTgxYTY0NDI5YWVhMmM2OTQ2MTdjNjI0YkBvNTUzNTU4LmluZ2VzdC5zZW50cnkuaW8vNTcyNjcwMQ=="
        ).decode("utf-8")
        # Capture exceptions on exit if onecodex CLI being invoked
        if os.path.basename(sys.argv[0]) in ["onecodex", "py.test"]:
            install_sys_hook = True
        else:
            install_sys_hook = False

        try:
            ipython_active = _setup_sentry_for_ipython()
            ignore_exceptions = []
            if ipython_active:
                # Patch the client to not raise too many exceptions in interactive environment
                # For now, we accept string and wildcard variables parsed from a special environment
                # variable. We can add support for hard-coded Exception classes here in the future as needed.
                ignore_exceptions = [
                    x
                    for x in os.environ.get("ONE_CODEX_SENTRY_IGNORE_EXCEPTIONS", "").split(",")
                    if x
                ]
            sentry_sdk.init(
                dsn=os.environ.get("ONE_CODEX_SENTRY_DSN", default_dsn),
                install_sys_hook=install_sys_hook,
                raise_send_errors=False,
                ignore_exceptions=ignore_exceptions,
                include_paths=[__name__.split(".", 1)[0]],
                release=__version__,
            )
            with sentry_sdk.configure_scope() as scope:
                if extra_context is None:
                    extra_context = {}
                if user_context is None:
                    user_context = {}
                scope.set_user(user_context)
                scope.set_extra("platform", platform.platform())
                if ipython_active:
                    scope.set_extra("ipython", True)

        except Exception:
            pass


def telemetry(fn):
    """Decorate CLI and other functions that need special Sentry client handling.

    This decorator is only required on functions that:
        1) Specifically catch and re-raise exceptions.
        2) Call sys.exit() directly.
    """

    @wraps(fn)
    def telemetry_wrapper(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except Exception as e:
            sentry_sdk.capture_exception(e)
            with sentry_sdk.configure_scope() as scope:
                scope.clear()
            raise e

    return telemetry_wrapper


def pretty_errors(fn):
    """Decorate CLI functions, catching errors and then calling sys.exit(1).

    This is intended for use with the CLI and scripts where we only use OneCodexException and its
    subclasses, e.g. ValidationError or UploadException.
    """

    @wraps(fn)
    def pretty_errors_wrapper(*args, **kwargs):
        try:
            fn(*args, **kwargs)
        except (OneCodexException, UploadException) as e:
            sys.stderr.write("\nERROR: {}".format(e))
            sys.stderr.write(
                "\nPlease feel free to contact us for help at support@onecodex.com\n\n"
            )
            sys.exit(1)

    return pretty_errors_wrapper


def snake_case(input_string):
    s1 = re.sub("(.)([A-Z][a-z]+)", r"\1_\2", input_string)
    return re.sub("([a-z0-9])([A-Z])", r"\1_\2", s1).lower()


def run_via_threadpool(fn, iterable, fn_kwargs, max_threads=1, graceful_exit=False):
    if max_threads == 1:
        for item in iterable:
            fn(item, **fn_kwargs)
        return

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = {executor.submit(fn, x, **fn_kwargs) for x in iterable}

        try:
            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                except Exception as e:
                    raise e
                else:
                    pass
        except KeyboardInterrupt as k:
            if not graceful_exit:
                executor._threads.clear()
                concurrent.futures.thread._threads_queues.clear()
            raise k


def progressbar(*args, **kwargs):
    bar = click.progressbar(*args, **kwargs)
    bar._update = bar.update

    def update(self, value):
        if getattr(self, "canceled", False) is True:
            self.label = "Canceling..."
            if self.eta_known:
                self.eta_known = 0
                self.pos = 0
                self._update(1)
        elif self.pct >= 0.98 and self.label.startswith("Uploading"):
            self.label = "Finalizing..."
            self._update(value)
        else:
            self._update(value)

    bar.update = partial(update, bar)
    return bar


# this lets us turn off the click progressbar context manager and is python2 compatible
# https://stackoverflow.com/questions/45187286/how-do-i-write-a-null-no-op-contextmanager-in-python
class FakeProgressBar(object):
    pct = 0
    label = ""

    def __init__(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def finish(self):
        pass

    def update(self, size):
        pass


def click_path_autocomplete_helper(ctx, args, incomplete, filename=True, directory=True):
    """Suggest paths to complete a partially typed filename or directory.

    Parameters
    ----------
    ctx : `click.Context`
    args : `list`
    incomplete : `str`
        A partially typed path to a file or directory.
    filename : `bool`
        If True, include filenames in the list of suggestions.
    directory : `bool`
        If True, include directories in the list of suggestions.

    Returns
    -------
    A list of suggestions for completing the path.
    """
    dir_name = os.path.dirname(incomplete)
    file_name = os.path.basename(incomplete)

    if dir_name.startswith("~"):
        abs_dir_name = os.path.expanduser(dir_name)
    else:
        abs_dir_name = os.path.abspath(dir_name)

    if not os.path.exists(abs_dir_name):
        return []

    suggestions = []

    for item in os.listdir(abs_dir_name):
        if item.startswith(file_name):
            abs_item_path = os.path.join(abs_dir_name, item)
            item_path = os.path.join(dir_name, item)

            if directory and os.path.isdir(abs_item_path):
                suggestions.append(item_path)
            elif filename and os.path.isfile(abs_item_path):
                suggestions.append(item_path)

    return suggestions


def is_continuous(series):
    import pandas as pd

    return not (
        pd.api.types.is_bool_dtype(series)
        or pd.api.types.is_categorical_dtype(series)  # noqa
        or pd.api.types.is_object_dtype(series)  # noqa
    )


def has_missing_values(dataframe_or_series):
    return dataframe_or_series.isnull().values.any()
