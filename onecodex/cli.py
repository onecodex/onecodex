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
import datetime

from onecodex.utils import (cli_resource_fetcher, download_file_helper,
                            valid_api_key, OPTION_HELP, pprint,
                            warn_if_insecure_platform, is_simplejson_installed,
                            warn_simplejson, telemetry)
from onecodex.api import Api
from onecodex.exceptions import ValidationWarning, ValidationError, UploadException
from onecodex.auth import _login, _logout, _remove_creds, _silent_login
from onecodex.scripts import filter_reads
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
@click.option("--telemetry/--no-telemetry", is_flag=True, default=True,
              help=OPTION_HELP['telemetry'])
@click.version_option(version=__version__)
@click.pass_context
@telemetry
def onecodex(ctx, api_key, no_pprint, verbose, telemetry):
    """One Codex v1 API command line interface"""
    # set up the context for sub commands
    click.Context.get_usage = click.Context.get_help
    ctx.obj = {}
    ctx.obj['API_KEY'] = api_key
    ctx.obj['NOPPRINT'] = no_pprint
    ctx.obj['TELEMETRY'] = telemetry

    if verbose:
        log.setLevel(logging.INFO)

    # Show a warning if simplejson is installed
    if is_simplejson_installed():
        warn_simplejson()

    # create the api
    no_api_subcommands = ["login", "logout"]
    if ctx.invoked_subcommand not in no_api_subcommands:
        if api_key is not None:
            ctx.obj['API'] = Api(cache_schema=True,
                                 api_key=api_key, telemetry=telemetry)
        else:
            # try and find it
            api_key = _silent_login()
            if api_key is not None:
                ctx.obj['API'] = Api(cache_schema=True, api_key=api_key, telemetry=telemetry)
            else:
                click.echo("No One Codex API credentials found - running anonymously", err=True)
                ctx.obj['API'] = Api(cache_schema=True, api_key='', telemetry=telemetry)

    # handle checking insecure platform, we let upload command do it by itself
    if ctx.invoked_subcommand != "upload":
        warn_if_insecure_platform()


@onecodex.group('scripts', help='Assorted utility scripts')
def scripts():
    pass


scripts.add_command(filter_reads.cli, 'filter_reads')


# resources
@onecodex.command('analyses')
@click.argument('analyses', nargs=-1, required=False)
@click.pass_context
@telemetry
def analyses(ctx, analyses):
    """Retrieve performed analyses"""
    cli_resource_fetcher(ctx, "analyses", analyses)


@onecodex.command('classifications')
@click.option("--read-level", 'readlevel', is_flag=True,
              help=OPTION_HELP['readlevel'])
@click.option("--read-level-path", 'readlevel_path', type=click.Path(),
              default="./", help=OPTION_HELP['readlevel_path'])
@click.option("--results", 'results', is_flag=True,
              help=OPTION_HELP['results'])
@click.pass_context
@click.argument('classifications', nargs=-1, required=False)
@telemetry
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
            if not classification:
                log.error('Could not find classification {} (404 status code)'
                          .format(classifications[0]))
                return
            results = classification.results(json=True)
            pprint(results, ctx.obj['NOPPRINT'])

    # fetch the readlevel
    elif readlevel is not None and not results:
        if len(classifications) != 1:
            log.error("Can only request read-level data on one Classification at a time")
        else:
            classification = ctx.obj['API'].Classifications.get(classifications[0])
            if not classification:
                log.error('Could not find classification {} (404 status code)'
                          .format(classifications[0]))
                return
            tsv_url = classification.readlevel()['url']
            log.info("Downloading tsv data from: {}".format(tsv_url))
            download_file_helper(tsv_url, readlevel_path)

    # both given -- complain
    else:
        log.error("Can only request one of read-level data or results data at a time")


@onecodex.command('panels')
@click.pass_context
@click.argument('panels', nargs=-1, required=False)
@telemetry
def panels(ctx, panels):
    """Retrieve performed in silico panel results"""
    cli_resource_fetcher(ctx, "panels", panels)


@onecodex.command('samples')
@click.pass_context
@click.argument('samples', nargs=-1, required=False)
@telemetry
def samples(ctx, samples):
    """Retrieve uploaded samples"""
    cli_resource_fetcher(ctx, "samples", samples)


# utilites
@onecodex.command('upload')
@click.option('--max-threads', default=4,
              help=OPTION_HELP['max_threads'], metavar='<int:threads>')
@click.argument('files', nargs=-1, required=False, type=click.Path(exists=True))
@click.option('--forward', type=click.Path(exists=True),
              help=OPTION_HELP['forward'])
@click.option('--reverse', type=click.Path(exists=True),
              help=OPTION_HELP['reverse'])
@click.option('--clean', is_flag=True, help=OPTION_HELP['clean'], default=False)
@click.option('--do-not-interleave', 'no_interleave', is_flag=True, help=OPTION_HELP['interleave'],
              default=False)
@click.option('--prompt/--no-prompt', is_flag=True, help=OPTION_HELP['prompt'], default=True)
@click.option('--validate/--do-not-validate', is_flag=True, help=OPTION_HELP['validate'],
              default=True)
@click.option('--tags', type=str)
@click.option("--metadata", '-md', multiple=True)
@click.pass_context
@telemetry
def upload(ctx, files, max_threads, clean, no_interleave, prompt, validate,
           forward, reverse, tags, metadata):
    """Upload a FASTA or FASTQ (optionally gzip'd) to One Codex"""

    tag_array = []
    if tags:
        tag_array = [tag.strip() for tag in tags.split(',')]

    metadata_opts = {}
    if metadata:
        for metadata_kv in metadata:
            split_metadata = metadata_kv.split('=')
            if len(split_metadata) > 1:
                metadata_value = '='.join(split_metadata[1:])
                metadata_opts[split_metadata[0].lower()] = metadata_value

    if (forward or reverse) and not (forward and reverse):
        click.echo('You must specify both forward and reverse files', err=True)
        sys.exit(1)
    if forward and reverse:
        if len(files) > 0:
            click.echo('You may not pass a FILES argument when using the '
                       ' --forward and --reverse options.', err=True)
            sys.exit(1)
        files = [(forward, reverse)]
        no_interleave = True
    if len(files) == 0:
        click.echo(ctx.get_help())
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
        sample_uuids = ctx.obj['API'].Samples.upload(files, threads=max_threads, validate=validate)
        update_sample_tags_and_metadata(sample_uuids, ctx, tag_array, metadata_opts)

    except ValidationWarning as e:
        sys.stderr.write('\nERROR: {}. {}'.format(
            e, 'Running with the --clean flag will suppress this error.'
        ))
        sys.exit(1)
    except (ValidationError, UploadException, Exception) as e:
        # TODO: Some day improve specific other exception error messages, e.g., gzip CRC IOError
        sys.stderr.write('\nERROR: {}'.format(e))
        sys.stderr.write('\nPlease feel free to contact us for help at help@onecodex.com\n\n')
        sys.exit(1)


def update_sample_tags_and_metadata(sample_uuids, ctx, tag_array, metadata_opts):
    if tag_array or metadata_opts:
        for uuid in sample_uuids:
            sample = ctx.obj['API'].Samples.get(uuid)
            if tag_array:
                new_tag_array = []
                for tag in tag_array:
                    existing_tags = ctx.obj['API'].Tags.where(name=tag)
                    if existing_tags:
                        unsaved_tag = existing_tags[0]
                        new_tag_array.append(unsaved_tag)
                    else:
                        unsaved_tag = ctx.obj['API'].Tags(name=tag, sample=sample)
                        unsaved_tag.save()

                if new_tag_array:
                    potential_tags = sample.tags
                    potential_tags.extend(new_tag_array)
                    sample.tags = potential_tags
                    sample.save()

            if metadata_opts:
                for k, v in metadata_opts.iteritems():
                    parse_metadata(sample, k, v)
                sample.metadata.save()


def parse_metadata(sample, metadata_key, metadata_value):
    if is_custom_metadata(metadata_key):
        sample.metadata.custom[metadata_key] = metadata_value
    else:
        set_metadata_value(sample.metadata, metadata_key, metadata_value)


def set_metadata_value(metadata, metadata_key, metadata_value):
    metadata_type = is_settable_metadata_keys()[metadata_key]
    if metadata_type == str:
        setattr(metadata, metadata_key, metadata_value)
    elif metadata_type == datetime.datetime:
        datetime_vals = metadata_value.split(',')
        datetime_ints = list(map(int, datetime_vals))
        setattr(metadata, metadata_key, datetime.datetime(*datetime_ints))
    elif metadata_type == bool:
        bool_val = metadata_value.lower() in ['true', '1', 't', 'y', 'yes']
        setattr(metadata, metadata_key, bool_val)
    elif metadata_type == 'not_supported':
        warning = "'{}' is not currently supported from the CLI".format(metadata_key)
        warnings.warn(warning)


def is_settable_metadata_keys():
    # TODO - Eventually, this shouldn't be hard-coded.
    return {
        'name': str,
        'description': str,
        'starred': bool,
        'date_collected': datetime.datetime,
        'date_sequenced': datetime.datetime,
        'sequencer': 'not_supported',
        'location': 'not_supported',
        'location_lon': 'not_supported',
        'location_lat': 'not_supported',
        'location_string': 'not_supported',
        'platform': 'not_supported'
    }


def is_custom_metadata(metadata_key):
    if is_settable_metadata_keys().get(metadata_key):
        return False
    else:
        return True


@onecodex.command('login')
@click.pass_context
@telemetry
def login(ctx):
    """Add an API key (saved in ~/.onecodex)"""
    base_url = os.environ.get("ONE_CODEX_API_BASE", "https://app.onecodex.com")
    if not ctx.obj['API_KEY']:
        _login(base_url)
    else:
        email = _login(base_url, api_key=ctx.obj['API_KEY'])
        ocx = Api(cache_schema=True, api_key=ctx.obj['API_KEY'], telemetry=ctx.obj['TELEMETRY'])

        # TODO: This should be protected or built in as a first class resource
        # with, e.g., connection error catching (it's not part of our formally documeted API at the moment)
        if ocx._client.Account.instances()['email'] != email:
            click.echo('Your login credentials do not match the provided email!', err=True)
            _remove_creds()
            sys.exit(1)


@onecodex.command('logout')
@click.pass_context
@telemetry
def logout(ctx):
    """Delete your API key (saved in ~/.onecodex)"""
    _logout()
