import base64
from click import BadParameter, Context, echo
from functools import wraps
import json
import logging
import os
import platform
import re
import requests
import sys

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse

from onecodex.vendored.potion_client.converter import PotionJSONEncoder
from onecodex.exceptions import OneCodexException, UploadException
from onecodex.version import __version__


log = logging.getLogger(__name__)
cli_log = logging.getLogger("onecodex.cli")

OPTION_HELP = {
    'api_key': 'Manually provide a One Codex API key',
    'no_pprint': 'Do not pretty-print JSON responses',
    'threads': 'Do not use multiple background threads to upload files',  # noqa
    'max_threads': 'Specify a different max # of upload threads (defaults to 4)',  # noqa
    'verbose': 'Log extra information to STDERR',
    'results': 'Get a JSON array of the metagenomic classification results table',
    'readlevel': 'Get the read-level data as a .tsv file',
    'readlevel_path': ('Output path or directory for a .tsv file with the raw '
                       'read-level results. Defaults to original filename in '
                       'the current working directory'),
    'clean': ("Automatically clean up FASTX records during upload. This removes tabs from "
              "headers and converts U to T. If this isn't passed, these will cause errors."),
    'interleave': ("Do not automatically interleave paired end files. Note this normally happens "
                   "during the upload, does not require additional disk space, and does not change "
                   "the files themselves."),
    'prompt': ("Manually prompt about automatic paired file interleaving. Setting --no-prompt "
               "will allow running without any user intervention, e.g. in a script."),
    'validate': ("Do not validate the FASTA/Q file before uploading. Incompatible with automatic "
                 "paired end interleaving (NOT RECOMMENDED)."),
    'telemetry': 'Send errors to One Codex?',
    'forward': 'Specify a forward reads file',
    'reverse': 'Specify a reverse reads file',
    'tag': ('NOTE: We recommend invoking the upload command separately for each '
            'sample to add sample-specific tags.\n\nAdd one or more tags to '
            'all uploaded samples e.g., `onecodex upload --tag "Cohort A" $FILE`.'),
    'metadata': ('NOTE: We recommend invoking the upload command separately for each '
                 'sample to add sample-specific metadata.\n\nAdd one or more '
                 'metadata attributes to all uploaded samples, '
                 'e.g. `onecodex upload --metadata starred=true --metadata '
                 'platform="Illumina MiSeq" $FILE`. '),
    'project': 'Provide the name, short name, or alphanumeric UUID of a ' \
               'project to automatically add the samples to that project on ' \
               'upload. NOTE: does not currently support adding a sample to ' \
               'a public project. Projects are searched by UUID, then name, ' \
               'then short name in that order.'
}

SUPPORTED_EXTENSIONS = ["fa", "fasta", "fq", "fastq",
                        "fa.gz", "fasta.gz", "fq.gz", "fastq.gz",
                        "fa.gzip", "fasta.gzip", "fq.gzip", "fastq.gzip"]


def valid_api_key(ctx, param, value):
    """
    Ensures an API has valid length (this is a click callback)
    """
    if value is not None and len(value) != 32:
        raise BadParameter(
            "API Key must be 32 characters long, not {}".format(str(len(value)))
        )
    else:
        return value


def pprint(j, no_pretty):
    """
    Prints as formatted JSON
    """
    if not no_pretty:
        echo(json.dumps(j, cls=PotionJSONEncoder, sort_keys=True,
                        indent=4, separators=(',', ': ')))
    else:
        echo(j)


def cli_resource_fetcher(ctx, resource, uris, print_results=True):
    try:
        # analyses is passed, want Analyses
        resource_name = resource[0].upper() + resource[1:]
        if len(uris) == 0:

            # if non given fetch all
            cli_log.debug("No %s IDs given, fetching all...", resource_name)
            instances = getattr(ctx.obj['API'], resource_name).all()
            cli_log.debug("Fetched %i %ss", len(instances), resource)
            objs_to_return = [x._resource._properties for x in instances]
        else:
            uris = list(set(uris))
            cli_log.debug("Fetching %s: %s", resource_name, ",".join(uris))

            instances = []
            for uri in uris:
                try:
                    instance = getattr(ctx.obj['API'], resource_name).get(uri)
                    if instance is not None:
                        instances.append(instance._resource._properties)
                    else:
                        cli_log.error('Could not find {} {} (404 status code)'.format(resource_name, uri))
                except requests.exceptions.HTTPError as e:
                    cli_log.error('Could not find %s %s (%d status code)'.format(
                        resource_name, uri, e.response.status_code
                    ))
            objs_to_return = instances

        if print_results:
            pprint(objs_to_return, ctx.obj['NOPPRINT'])
        else:
            return objs_to_return
    except requests.exceptions.HTTPError:
        echo("Failed to authenticate. Please check your API key "
             "or trying logging out and back in with `onecodex logout` "
             "and `onecodex login`.")


def is_insecure_platform():
    """
    Checks if the current system is missing an SSLContext object
    """
    v = sys.version_info
    if v.major == 3:
        return False  # Python 2 issue

    if v.major == 2 and v.minor == 7 and v.micro >= 9:
        return False  # >= 2.7.9 includes the new SSL updates

    try:
        import OpenSSL  # noqa
        import ndg  # noqa
        import pyasn1  # noqa
    except ImportError:
        pass
    return True


def warn_if_insecure_platform():
    """
    Produces a nice message if SSLContext object is not available.
    Also returns True -> platform is insecure
                 False -> platform is secure
    """
    m = ("\n"
         "######################################################################################\n"  # noqa
         "#                                                                                    #\n"  # noqa
         "#  Your version of Python appears to be out of date and lack important security      #\n"  # noqa
         "#  features. Please update to Python >= 2.7.9 or `pip install requests[security]`.   #\n"  # noqa
         "#                                                                                    #\n"  # noqa
         "#  InsecurePlatformWarning: A true SSLContext object is not available. This          #\n"  # noqa
         "#  prevents urllib3 from configuring SSL appropriately and may cause certain         #\n"  # noqa
         "#  SSL connections to fail. For more information, see                                #\n"  # noqa
         "#  https://urllib3.readthedocs.org/en/latest/security.html#insecureplatformwarning.  #\n"  # noqa
         "#                                                                                    #\n"  # noqa
         "######################################################################################\n")  # noqa
    if is_insecure_platform():
        echo(m, err=True)
        return True
    else:
        cli_log.debug("Python SSLContext passed")
        return False


def get_download_dest(input_path, url):
    original_filename = urlparse(url).path.split("/")[-1]
    if os.path.isdir(input_path):
        local_full_path = os.path.join(input_path, original_filename)
    else:
        local_full_path = input_path
    return local_full_path


def download_file_helper(url, input_path):
    """
    Manages the chunked downloading of a file given an url
    """
    r = requests.get(url, stream=True)
    if r.status_code != 200:
        cli_log.error("Failed to download file: %s" % r.json()["message"])
    local_full_path = get_download_dest(input_path, r.url)
    original_filename = os.path.split(local_full_path)[-1]
    with open(local_full_path, 'wb') as f:
        echo("Downloading {}".format(original_filename), err=True)
        for chunk in r.iter_content(chunk_size=1024):
            if chunk:  # filter out keep-alive new chunks
                f.write(chunk)
                f.flush()
    pprint("Successfully downloaded %s to %s" % (original_filename, local_full_path),  # noqa
           True)


def check_for_allowed_file(f):
    """
    Checks a file extension against a list of seq file exts
    """
    for ext in SUPPORTED_EXTENSIONS:
        if f.endswith(ext):
            return True
    log.error("Failed upload: Not an allowed file extension: %s", f)
    raise SystemExit


def collapse_user(fp):
    """
    Converts a path back to ~/ from expanduser()
    """
    home_dir = os.path.expanduser("~")
    abs_path = os.path.abspath(fp)
    return abs_path.replace(home_dir, "~")


def _setup_sentry_for_ipython(client):
    from IPython import get_ipython
    ip = get_ipython()

    def custom_exc(shell, etype, evalue, tb, tb_offset=None):
        # Show original error
        shell.showtraceback((etype, evalue, tb), tb_offset=tb_offset)

        # Then send to Sentry
        client.captureException()

    if ip is not None:
        ip.set_custom_exc((Exception,), custom_exc)

    # Finally, patch the client to not raise too many exceptions in interactive environment
    # For now, we accept string and wildcard variables parsed from a special environment
    # variable. We can add support for hard-coded Exception classes here in the future as needed.
    client.ignore_exceptions = [
        x for x in os.environ.get('ONE_CODEX_SENTRY_IGNORE_EXCEPTIONS', '').split(',') if x
    ]


def get_raven_client(user_context=None, extra_context=None):
    if os.environ.get('ONE_CODEX_NO_TELEMETRY') is None:
        key = base64.b64decode(
            b'NmFlMjMwYWY4NjI5NDg3NmEyYzYwYjZjNDhhZDJiYzI6ZTMyZmYwZTVhNjUwNGQ5NGJhODc0NWZlMmU1ZjNmZjA='
        ).decode('utf-8')

        # Set Client params
        # Capture exceptions on exit if onecodex CLI being invoked
        if os.path.basename(sys.argv[0]) in ['onecodex', 'py.test']:
            install_sys_hook = True
        else:
            install_sys_hook = False

        try:
            from raven import Client as RavenClient
            client = RavenClient(
                dsn=os.environ.get('ONE_CODEX_SENTRY_DSN',
                                   'https://{}@sentry.onecodex.com/9'.format(key)),
                install_sys_hook=install_sys_hook,
                raise_send_errors=False,
                ignore_exceptions=[],
                include_paths=[__name__.split('.', 1)[0]],
                release=__version__
            )

            if extra_context is None:
                extra_context = {}
            if user_context is None:
                user_context = {}

            try:
                _setup_sentry_for_ipython(client)
                extra_context['ipython'] = True
            except Exception:
                pass

            extra_context['platform'] = platform.platform()
            client.user_context(user_context)
            client.extra_context(extra_context)
            return client

        except Exception:
            return


def telemetry(fn):
    """
    Decorator for CLI and other functions that need special Sentry client handling.
    This function is only required for functions that may exit *before* we set up
    the ._raven_client object on the Api instance *or* that specifically catch and re-raise
    exceptions or call sys.exit directly.

    Note that this also overwrites verbose Raven logs on exit ("Sentry is waiting to send..."),
    see https://github.com/getsentry/raven-python/issues/904 for more details.
    """

    @wraps(fn)
    def telemetry_wrapper(*args, **kwargs):
        # By default, do not instantiate a client,
        # and inherit the telemetry settings passed
        client = None
        if len(args) > 0 and isinstance(args[0], Context):
            ctx = args[0]

            # First try to get off API obj
            if ctx.obj and ctx.obj.get('API') is not None:
                client = ctx.obj['API']._raven_client

            # Else try to see if telemetry param is set
            elif ctx.params.get('telemetry', False):
                client = get_raven_client()

            # Finally check for the ctx.obj['TELEMETRY'] bool
            elif ctx.obj and ctx.obj.get('TELEMETRY', False):
                client = get_raven_client()

        try:
            return fn(*args, **kwargs)
        except Exception as e:
            if client:
                client.captureException()
                client.context.clear()
                sys.stdout = StringIO()
            raise e

    return telemetry_wrapper


def pretty_errors(fn):
    """
    Decorator for the CLI for catching errors and then calling sys.exit(1).

    For now, this is intended for use with the CLI and scripts where we only use
    OneCodexExceptions (incl. ValidationError) and UploadException.
    """
    @wraps(fn)
    def pretty_errors_wrapper(*args, **kwargs):
        try:
            fn(*args, **kwargs)
        except (OneCodexException, UploadException) as e:
            sys.stderr.write('\nERROR: {}'.format(e))
            sys.stderr.write('\nPlease feel free to contact us for help at help@onecodex.com\n\n')
            sys.exit(1)
    return pretty_errors_wrapper


def snake_case(input_string):
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', input_string)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()
