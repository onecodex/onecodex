"""
utils.py
author: @mbiokyle29
"""
import json
import logging
import os
import sys
try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse


import requests
from click import BadParameter, echo
from potion_client.converter import PotionJSONEncoder

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
}

SUPPORTED_EXTENSIONS = ["fa", "fasta", "fq", "fastq",
                        "fa.gz", "fasta.gz", "fq.gz", "fastq.gz",
                        "fa.gzip", "fasta.gzip", "fq.gzip", "fastq.gzip"]


def valid_api_key(ctx, param, value):
    """
    Ensures an API has valid length (this is a click callback)"""
    if value is not None and len(value) != 32:
        l = len(value)
        raise BadParameter(
            "API Key must be 32 characters long, not {}".format(str(l)))
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


def cli_resource_fetcher(ctx, resource, uris):
    """Helper method to parse CLI args in API calls
    """
    try:
        _cli_resource_fetcher(ctx, resource, uris)
    except requests.exceptions.HTTPError:
        echo("Failed to authenticate. Please check your API key "
             "or trying logging out and back in with `onecodex logout` "
             "and `onecodex login`.")


def _cli_resource_fetcher(ctx, resource, uris):
    # analyses is passed, want Analyses
    resource_name = resource[0].upper() + resource[1:]
    if len(uris) == 0:

        # if non given fetch all
        cli_log.info("No %s IDs given, fetching all...", resource_name)
        instances = getattr(ctx.obj['API'], resource_name).all()
        cli_log.info("Fetched %i %ss", len(instances), resource)
        pprint([x._resource._properties for x in instances], ctx.obj['NOPPRINT'])
    else:
        uris = list(set(uris))
        cli_log.info("Fetching %s: %s", resource_name, ",".join(uris))

        instances = []
        for uri in uris:
            try:
                instance = getattr(ctx.obj['API'], resource_name).get(uri)
                instances.append(instance._resource._properties)
            except requests.exceptions.HTTPError as e:
                cli_log.error("Could not find %s %s (%d status code)",
                              resource_name, uri, e.response.status_code)
        pprint(instances, ctx.obj['NOPPRINT'])


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
        cli_log.info("Python SSLContext passed")
        return False


def download_file_helper(url, input_path):
    """
    Manages the chunked downloading of a file given an url
    """
    r = requests.get(url, stream=True)
    if r.status_code != 200:
        cli_log.error("Failed to download file: %s" % r.json()["message"])
    original_filename = urlparse(r.url).path.split("/")[-1]
    if os.path.isdir(input_path):
        local_full_path = os.path.join(input_path, original_filename)
    else:
        local_full_path = input_path
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
