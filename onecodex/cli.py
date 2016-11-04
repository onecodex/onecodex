#!/usr/bin/env python
"""
cli.py
author: @mbiokyle29
"""
from __future__ import print_function
import os
import logging

import click

from onecodex.utils import (cli_resource_fetcher, download_file_helper,
                            valid_api_key, OPTION_HELP, pprint,
                            warn_if_insecure_platform)
from onecodex.api import Api
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
                click.echo("No One Codex API key available - running anonymously", err=True)
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
@click.option("--raw", 'raw', is_flag=True,
              help=OPTION_HELP['raw'])
@click.option("--raw-path", 'raw_path',
              default="./", help=OPTION_HELP['raw_path'])
@click.option("--results", 'results', is_flag=True,
              help=OPTION_HELP['results'])
@click.pass_context
@click.argument('classifications', nargs=-1, required=False)
def classifications(ctx, classifications, results, raw, raw_path):
    """Retrieve performed metagenomic classifications"""

    # basic operation -- just print
    if not raw and not results:
        cli_resource_fetcher(ctx, "classifications", classifications)

    # fetch the results
    elif not raw and results:
        if len(classifications) != 1:
            log.error("Can only request results data on one Classification at a time")
        else:
            classification = ctx.obj['API'].Classifications.get(classifications[0])
            results = classification.results(json=True)
            pprint(results, ctx.obj['NOPPRINT'])

    # fetch the raw
    elif raw is not None and not results:

        if len(classifications) != 1:
            log.error("Can only request raw data on one Classification at a time")

        else:
            classification = ctx.obj['API'].Classifications.get(classifications[0])
            tsv_url = classification.readlevel()['url']
            log.info("Downloading tsv data from: {}".format(tsv_url))
            download_file_helper(tsv_url, raw_path)

    # both given -- complain
    else:
        log.error("Can only request one of raw data or results data at a time")


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
@click.option("--no-threads", is_flag=True,
              help=OPTION_HELP['threads'],
              default=False)
@click.option("--max-threads", default=4,
              help=OPTION_HELP['max_threads'], metavar="<int:threads>")
@click.argument('files', nargs=-1, required=False,
                type=click.Path(exists=True))
@click.pass_context
def upload(ctx, files, no_threads, max_threads):
    """Upload a FASTA or FASTQ (optionally gzip'd) to One Codex"""
    if not no_threads:
        ctx.obj['API'].Samples.upload(files, threads=max_threads)
    else:
        ctx.obj['API'].Samples.upload(files, threads=1)


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
