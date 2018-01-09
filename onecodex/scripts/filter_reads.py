import click
import csv
import gzip
import os

from onecodex.exceptions import OneCodexException, ValidationError
from onecodex.lib.inline_validator import FASTXTranslator
from onecodex.utils import download_file_helper, get_download_dest, pretty_errors


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


def write_fastx_record(record, handler):
    if len(record) == 2:
        record_str = '>{}\n{}'
    elif len(record) == 4:
        record_str = '@{}\n{}\n{}\n{}'
    else:
        raise OneCodexException('Unknown FASTX record format', record)
    handler.write(record_str.format(*record))


@click.command(help='Filter a FASTX file based on the taxonomic results from a CLASSIFICATION_ID')
@click.argument('classification_id')
@click.argument('fastx', type=click.Path())
@click.option('-t', '--tax-id', required=True, multiple=True,
              help='Filter to reads mapping to tax IDs. May be passed multiple times.')
@click.option('-r', '--reverse', type=click.Path(), help='The reverse (R2) '
              'read file, optionally')
@click.option('--split-pairs/--keep-pairs', default=False, help='Keep only '
              'the read pair member that matches the list of tax ID\'s')
@click.option('-o', '--out', default='.', type=click.Path(), help='Where '
              'to put the filtered outputs')
@click.pass_context
@pretty_errors
def cli(ctx, classification_id, fastx, reverse, tax_id, split_pairs, out):
    tax_ids = tax_id  # rename
    if not len(tax_ids):
        raise OneCodexException('You must supply at least one tax ID')

    classification = ctx.obj['API'].Classifications.get(classification_id)
    if classification is None:
        raise ValidationError('Classification {} not found.'.format(classification_id))

    tsv_url = classification.readlevel()['url']
    readlevel_path = get_download_dest('./', tsv_url)
    if not os.path.exists(readlevel_path):
        download_file_helper(tsv_url, './')
    else:
        click.echo('Using cached read-level results: {}'
                   .format(readlevel_path), err=True)

    filtered_rows = []
    tsv_row_count = 0
    with gzip.open(readlevel_path, 'rt') as tsv:
        try:
            tsv_row_count = get_record_count(tsv) - 1  # discount header line
        except EOFError:
            click.echo('\nWe encountered an error while processing the read '
                       'level results. Please delete {} and try again.'
                       .format(readlevel_path), err=True)
            raise
        else:
            tsv.seek(0)
            reader = csv.DictReader(tsv, delimiter='\t')
            click.echo('Selecting results matching tax ID(s): {}'
                       .format(', '.join(tax_ids)), err=True)
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

    save_msg = 'Saving filtered reads: {}'.format(filtered_filename)
    if reverse:
        save_msg += ' and {}'.format(rev_filtered_filename)
    click.echo(save_msg, err=True)

    counter = 0
    if reverse:
        with open(fastx, 'rb') as fastx_file, \
                open(reverse, 'rb') as reverse_file, \
                open(filtered_filename, 'wb') as out_file, \
                open(rev_filtered_filename, 'wb') as rev_out_file:
            if split_pairs:
                for fwd, rev in FASTXTranslator(fastx_file, reverse_file,
                                                validate=False):
                    if counter in filtered_rows:
                        out_file.write(fwd)
                    if (counter + 1) in filtered_rows:
                        rev_out_file.write(rev)
                    counter += 2
            else:
                for fwd, rev in FASTXTranslator(fastx_file, reverse_file,
                                                validate=False):
                    if counter in filtered_rows or \
                       (counter + 1) in filtered_rows:
                        out_file.write(fwd)
                        rev_out_file.write(rev)
                    counter += 2
    else:
        with open(fastx, 'rb') as fastx_file, \
                open(filtered_filename, 'wb') as out_file:
            for seq in FASTXTranslator(fastx_file, validate=False):
                if counter in filtered_rows:
                    out_file.write(seq)
                counter += 1
