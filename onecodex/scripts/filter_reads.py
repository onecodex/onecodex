import click
import csv
import gzip
import os

from onecodex.exceptions import OneCodexException, ValidationError
from onecodex.lib.inline_validator import FASTXTranslator
from onecodex.utils import download_file_helper, get_download_dest


def with_progress_bar(length, ix, *args, **kwargs):
    label = kwargs.pop('label', None)
    bar_kwargs = {
        'length': length
    }
    if label:
        bar_kwargs['label'] = label
    with click.progressbar(**bar_kwargs) as bar:
        return ix(cb=bar.update, *args, **kwargs)


def filter_rows_by_taxid(results, tax_ids, cb=None):
    filtered_rows = []
    for i, row in enumerate(results):
        if row['Tax ID'] in tax_ids:
            filtered_rows.append(i)
        if i % 1000 == 0 and cb:
            cb(1000)
    return filtered_rows


def get_record_count(fp):
    counter = 0
    for _ in fp:
        counter += 1
    return counter


def get_filtered_filename(filepath):
    filename = os.path.split(filepath)[-1]
    prefix, ext = os.path.splitext(filename.rstrip('.gz')
                                           .rstrip('gzip'))
    return '{}.filtered{}'.format(prefix, ext), ext


def fetch_readlevel_results(api, classification_id):
    classification = api.Classifications.get(classification_id)
    if not classification:
        raise ValidationError('Classification {} '
                              'not found'.format(classification_id))
    tsv_url = classification.readlevel()['url']
    download_file_helper(tsv_url, './')
    return get_download_dest('./', tsv_url)


def write_fastx_record(record, handler):
    if len(record) == 2:
        record_str = '>{}\n{}'
    elif len(record) == 4:
        record_str = '@{}\n{}\n{}\n{}'
    else:
        raise OneCodexException('Unknown FASTX record format', record)
    handler.write(record_str.format(*record))


@click.command(help='Filters read files based on tax ID\'s')
@click.argument('classification')
@click.argument('fastx', type=click.Path())
@click.option('-t', '--tax-ids', required=True, help='A comma-delimited list '
              'of tax ID\'s to retain')
@click.option('-r', '--reverse', type=click.Path(), help='The reverse (R2) '
              'read file, optionally')
@click.option('--split-pairs/--keep-pairs', default=False, help='Keep only '
              'the read pair member that matches the list of tax ID\'s')
@click.option('-o', '--out', default='.', type=click.Path(), help='Where '
              'to put the filtered outputs')
@click.pass_context
def cli(ctx, classification, fastx, reverse, tax_ids, split_pairs, out):
    tax_ids = tax_ids.split(',')
    if not len(tax_ids):
        raise OneCodexException('You must supply at least one tax ID')

    readlevel_path = fastx \
        .replace('_R1', '') \
        .replace('_R2', '') \
        .rstrip('.gz') \
        .rstrip('.gzip') + '.gz.results.tsv.gz'
    if not os.path.exists(readlevel_path):
        click.echo('Downloading read-level results')
        readlevel_path = fetch_readlevel_results(ctx.obj['API'],
                                                 classification)
    else:
        click.echo('Using cached read-level results from {}'
                   .format(readlevel_path))

    filtered_rows = []
    tsv_row_count = 0
    with gzip.open(readlevel_path, 'rt') as tsv:
        try:
            click.echo('Getting result count')
            tsv_row_count = get_record_count(tsv) - 1  # discount header line
        except EOFError:
            click.echo('\nWe encountered an error while processing the read '
                       'level results. Please delete {} and try again.'
                       .format(readlevel_path), err=True)
            raise
        else:
            tsv.seek(0)
            reader = csv.DictReader(tsv, delimiter='\t')
            click.echo('Selecting results matching tax ID(s) {}'
                       .format(','.join(tax_ids)))
            filtered_rows = with_progress_bar(
                tsv_row_count,
                filter_rows_by_taxid,
                reader,
                tax_ids
            )

    filtered_filename = get_filtered_filename(fastx)[0]
    filtered_filename = os.path.join(out, filtered_filename)
    if reverse:
        rev_filtered_filename = get_filtered_filename(reverse)[0]
        rev_filtered_filename = os.path.join(out, rev_filtered_filename)

    click.echo('Getting FASTX record count')
    fastx_record_count = 0
    with open(fastx, 'rb') as fastx_file:
        try:
            fastx_record_count = get_record_count(
                FASTXTranslator(fastx_file, validate=False))
        except ValidationError as e:
            raise OneCodexException(e.message)

    if reverse:
        fastx_record_count = fastx_record_count * 2

    if tsv_row_count != fastx_record_count:
        os.remove(readlevel_path)
        raise OneCodexException('The supplied file has a different number of '
                                'records than the requested Classification')

    if reverse:
        click.echo('Filtering {} and {}'.format(fastx, reverse))
        with open(fastx, 'rb') as fastx_file, \
                open(reverse, 'rb') as reverse_file, \
                open(filtered_filename, 'wb') as out_file, \
                open(rev_filtered_filename, 'wb') as rev_out_file:
            if split_pairs:
                counter = 0
                for fwd, rev in FASTXTranslator(fastx_file, reverse_file,
                                                validate=False):
                    if counter in filtered_rows:
                        out_file.write(fwd)
                    if (counter + 1) in filtered_rows:
                        rev_out_file.write(rev)
                    counter += 2
            else:
                counter = 0
                for fwd, rev in FASTXTranslator(fastx_file, reverse_file,
                                                validate=False):
                    if counter in filtered_rows or \
                       (counter + 1) in filtered_rows:
                        out_file.write(fwd)
                        rev_out_file.write(rev)
                    counter += 2
    else:
        click.echo('Filtering {}'.format(fastx))
        with open(fastx, 'rb') as fastx_file, \
                open(filtered_filename, 'wb') as out_file:
            counter = 0
            for seq in FASTXTranslator(fastx_file, validate=False):
                if counter in filtered_rows:
                    out_file.write(seq)
                counter += 1
