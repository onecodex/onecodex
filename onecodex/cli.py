#!/usr/bin/env python
"""
cli.py
author: @mbiokyle29
"""
from __future__ import print_function
import logging
import os
import re
import sys
import warnings

import click

from onecodex.utils import (cli_resource_fetcher, download_file_helper,
                            valid_api_key, OPTION_HELP, pprint,
                            warn_if_insecure_platform)
from onecodex.api import Api
from onecodex.exceptions import ValidationWarning, ValidationError, UploadException
from onecodex.auth import _login, _logout, _silent_login
from onecodex.version import __version__

# set the context for getting -h also
CONTEXT_SETTINGS = dict(
    help_option_names=['-h', '--help'],
)

# logging
log = logging.getLogger(__name__)
log.setLevel(logging.WARN)
log_formatter = logging.Formatter('One Codex CLI (%(levelname)s): %(message)s')
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.INFO)
stream_handler.setFormatter(log_formatter)
log.addHandler(stream_handler)


# options
@click.group(context_settings=CONTEXT_SETTINGS)
@click.option("--api-key", callback=valid_api_key,
              help=OPTION_HELP['api_key'], metavar="<str:api_key>")
@click.option("--no-pretty-print", 'no_pprint', is_flag=True,
              help=OPTION_HELP['no_pprint'])
@click.option("--verbose", "-v", is_flag=True,
              help=OPTION_HELP['verbose'])
@click.version_option(version=__version__)
@click.pass_context
def onecodex(ctx, api_key, no_pprint, verbose):
    """One Codex v1 API command line interface"""

    # set up the context for sub commands
    click.Context.get_usage = click.Context.get_help
    ctx.obj = {}
    ctx.obj['API_KEY'] = api_key
    ctx.obj['NOPPRINT'] = no_pprint

    if verbose:
        log.setLevel(logging.INFO)

    if os.environ.get("ONE_CODEX_API_BASE") is not None:
        base_url = os.environ.get("ONE_CODEX_API_BASE")
        click.echo("ALL REQUESTS GOING THROUGH: %s" % base_url, err=True)
    else:
        base_url = "https://app.onecodex.com"

    ctx.obj['BASE_URL'] = base_url

    # create the api
    no_api_subcommands = ["login", "logout"]
    if ctx.invoked_subcommand not in no_api_subcommands:
        if api_key is not None:
            ctx.obj['API'] = Api(base_url=base_url, extensions=False,
                                 cache_schema=True,
                                 api_key=api_key)
        else:
            # try and find it
            api_key = _silent_login()
            if api_key is not None:
                ctx.obj['API'] = Api(base_url=base_url, extensions=False,
                                     cache_schema=True, api_key=api_key)
            else:
                click.echo("No One Codex API key is available - running anonymously", err=True)
                ctx.obj['API'] = Api(base_url=base_url, extensions=False, cache_schema=True)

    # handle checking insecure platform, we let upload command do it by itself
    if ctx.invoked_subcommand != "upload":
        warn_if_insecure_platform()


# resources
@onecodex.command('analyses')
@click.argument('analyses', nargs=-1, required=False)
@click.pass_context
def analyses(ctx, analyses):
    """Retrieve performed analyses"""
    cli_resource_fetcher(ctx, "analyses", analyses)


@onecodex.command('classifications')
@click.option("--read-level", 'readlevel', is_flag=True,
              help=OPTION_HELP['readlevel'])
@click.option("--read-level-path", 'readlevel_path',
              default="./", help=OPTION_HELP['readlevel_path'])
@click.option("--results", 'results', is_flag=True,
              help=OPTION_HELP['results'])
@click.pass_context
@click.argument('classifications', nargs=-1, required=False)
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
            classification = ctx.obj['API'].Classifications.get(classifications[0])
            results = classification.results(json=True)
            pprint(results, ctx.obj['NOPPRINT'])

    # fetch the readlevel
    elif readlevel is not None and not results:
        if len(classifications) != 1:
            log.error("Can only request read-level data on one Classification at a time")

        else:
            classification = ctx.obj['API'].Classifications.get(classifications[0])
            tsv_url = classification.readlevel()['url']
            log.info("Downloading tsv data from: {}".format(tsv_url))
            download_file_helper(tsv_url, readlevel_path)

    # both given -- complain
    else:
        log.error("Can only request one of read-level data or results data at a time")


@onecodex.command('panels')
@click.pass_context
@click.argument('panels', nargs=-1, required=False)
def panels(ctx, panels):
    """Retrieve performed in silico panel results"""
    cli_resource_fetcher(ctx, "panels", panels)


@onecodex.command('samples')
@click.pass_context
@click.argument('samples', nargs=-1, required=False)
def samples(ctx, samples):
    """Retrieve uploaded samples"""
    cli_resource_fetcher(ctx, "samples", samples)


# utilites
@onecodex.command('upload')
@click.option('--max-threads', default=4,
              help=OPTION_HELP['max_threads'], metavar='<int:threads>')
@click.argument('files', nargs=-1, required=False, type=click.Path(exists=True))
@click.option('--clean', is_flag=True, help=OPTION_HELP['clean'], default=False)
@click.option('--do-not-interleave', 'no_interleave', is_flag=True, help=OPTION_HELP['interleave'],
              default=False)
@click.option('--prompt/--no-prompt', is_flag=True, help=OPTION_HELP['prompt'], default=True)
@click.pass_context
def upload(ctx, files, max_threads, clean, no_interleave, prompt):
    """Upload a FASTA or FASTQ (optionally gzip'd) to One Codex"""
    if len(files) == 0:
        print(ctx.get_help())
        return
    else:
        files = list(files)

    if not no_interleave:
        # "intelligently" find paired files and tuple them
        paired_files = []
        single_files = set(files)
        for filename in files:
            # convert "read 1" filenames into "read 2" and check that they exist; if they do
            # upload the files as a pair, autointerleaving them
            pair = re.sub('[._][Rr]1[._]', lambda x: x.group().replace('1', '2'), filename)
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
            pair_list = ''
            for p in paired_files:
                pair_list += '\n  {}  &  {}'.format(os.path.basename(p[0]), os.path.basename(p[1]))

            answer = click.confirm(
                'It appears there are paired files:{}\nInterleave them after upload?'.format(
                    pair_list
                ),
                default='Y'
            )
            if not answer:
                auto_pair = False

        if auto_pair:
            files = paired_files + list(single_files)

    if not clean:
        warnings.filterwarnings('error', category=ValidationWarning)

    try:
        # do the uploading
        ctx.obj['API'].Samples.upload(files, threads=max_threads)
    except ValidationWarning as e:
        sys.stderr.write('\nERROR: {}. {}'.format(
            e, 'Running with the --clean flag will suppress this error.'
        ))
        sys.exit(1)
    except (ValidationError, UploadException, Exception) as e:
        # TODO: Some day improve specific other exception error messages, e.g., gzip CRC IOError
        sys.stderr.write('\nERROR: {}'.format(e))
        sys.stderr.write('\nPlease feel free to contact us for help at help@onecodex.com')
        sys.exit(1)


@onecodex.command('login')
@click.pass_context
def login(ctx):
    """Add an API key (saved in ~/.onecodex)"""
    _login(ctx.obj['BASE_URL'], check_for_update=False)


@onecodex.command('logout')
@click.pass_context
def logout(ctx):
    """Delete your API key (saved in ~/.onecodex)"""
    _logout()
