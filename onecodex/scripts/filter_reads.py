import click
import csv
import gzip
import io
import os

from onecodex.auth import login_required
from onecodex.exceptions import OneCodexException
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


def make_taxonomy_dict(classification, parent=False):
    """
    Takes a classification data frame returned by the API and parses it
    into a dictionary mapping a tax_id to its children (or parent).
    Restricted to tax_id's that are represented in the classification
    results.
    """

    tax_id_map = {}

    if parent:
        for row in classification.results()['table']:
            if row['parent_tax_id'] is not None and \
               row['tax_id'] is not None:
                tax_id_map[row['tax_id']] = row['parent_tax_id']
    else:
        for row in classification.results()['table']:
            if row['parent_tax_id'] is not None and \
               row['tax_id'] is not None:
                try:
                    tax_id_map[row['parent_tax_id']].add(row['tax_id'])
                except KeyError:
                    tax_id_map[row['parent_tax_id']] = set([row['tax_id']])

    return tax_id_map


def recurse_taxonomy_map(tax_id_map, tax_id, parent=False):
    """
    Takes the output dict from make_taxonomy_map and an input tax_id
    and recurses either up or down through the tree to get /all/ children
    (or parents) of the given tax_id.
    """

    if parent:
        pass
    else:
        def _child_recurse(tax_id, visited):
            try:
                children = [tax_id] + list(tax_id_map[tax_id])
            except KeyError:
                children = [tax_id]

            for child in children:
                if child not in visited:
                    visited.append(child)
                    children.extend(_child_recurse(child, visited))

            return children

        return list(set(_child_recurse(tax_id, [])))


@click.command('filter_reads', help='Filter a FASTX file based on the taxonomic results from a CLASSIFICATION_ID')
@click.argument('classification_id')
@click.argument('fastx', type=click.Path())
@click.option('-t', '--tax-id', required=True, multiple=True,
              help='Filter reads mapping to tax IDs. May be passed multiple times.')
@click.option('--with-children', default=False, is_flag=True,
              help='Keep reads of child taxa, too. For example, all strains of E. coli')
@click.option('-r', '--reverse', type=click.Path(), help='The reverse (R2) '
              'read file, optionally')
@click.option('--split-pairs/--keep-pairs', default=False, help='Keep only '
              'the read pair member that matches the list of tax ID\'s')
@click.option('--exclude-reads', default=False, is_flag=True,
              help='Rather than keep reads matching reads, choosing this option will exclude them.')
@click.option('-o', '--out', default='.', type=click.Path(), help='Where '
              'to put the filtered outputs')
@click.pass_context
@pretty_errors
@login_required
def cli(ctx, classification_id, fastx, reverse, tax_id, with_children,
        split_pairs, exclude_reads, out):
    import skbio

    if not len(tax_id):
        raise OneCodexException('You must supply at least one tax ID')

    classification = ctx.obj['API'].Classifications.get(classification_id)
    if classification is None:
        raise OneCodexException('Classification {} not found.'.format(classification_id))

    if with_children:
        tax_id_map = make_taxonomy_dict(classification)

        new_tax_ids = []

        for t_id in tax_id:
            new_tax_ids.extend(recurse_taxonomy_map(tax_id_map, t_id))

        tax_id = list(set(new_tax_ids))

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
                       .format(', '.join(tax_id)), err=True)
            filtered_rows = with_progress_bar(
                tsv_row_count,
                filter_rows_by_taxid,
                reader,
                tax_id
            )

    filtered_filename, ext = get_filtered_filename(fastx)
    filtered_filename = os.path.join(out, filtered_filename)
    if reverse:
        rev_filtered_filename = get_filtered_filename(reverse)[0]
        rev_filtered_filename = os.path.join(out, rev_filtered_filename)

    if ext in {'.fa', '.fna', '.fasta'}:
        io_kwargs = {'format': 'fasta'}
    elif ext in {'.fq', '.fastq'}:
        io_kwargs = {'format': 'fastq', 'variant': 'illumina1.8'}
    else:
        raise OneCodexException(
            '{}: extension must be one of .fa, .fna, .fasta, .fq, .fastq'.format(fastx)
        )

    fastx_record_count = get_record_count(skbio.io.read(fastx, **io_kwargs))

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

    # skbio doesn't support built-in open() method in python2. must use io.open()
    counter = 0
    if reverse:
        with io.open(filtered_filename, 'w') as out_file, \
             io.open(rev_filtered_filename, 'w') as rev_out_file:  # noqa
            if split_pairs:
                for fwd, rev in zip(skbio.io.read(fastx, **io_kwargs),
                                    skbio.io.read(reverse, **io_kwargs)):
                    if exclude_reads:
                        if counter not in filtered_rows:
                            fwd.write(out_file, **io_kwargs)
                        if (counter + 1) not in filtered_rows:
                            rev.write(rev_out_file, **io_kwargs)
                    else:
                        if counter in filtered_rows:
                            fwd.write(out_file, **io_kwargs)
                        if (counter + 1) in filtered_rows:
                            rev.write(rev_out_file, **io_kwargs)
                    counter += 2
            else:
                for fwd, rev in zip(skbio.io.read(fastx, **io_kwargs),
                                    skbio.io.read(reverse, **io_kwargs)):
                    if exclude_reads:
                        if counter not in filtered_rows and \
                           (counter + 1) not in filtered_rows:
                            fwd.write(out_file, **io_kwargs)
                            rev.write(rev_out_file, **io_kwargs)
                    else:
                        if counter in filtered_rows or \
                           (counter + 1) in filtered_rows:
                            fwd.write(out_file, **io_kwargs)
                            rev.write(rev_out_file, **io_kwargs)
                    counter += 2
    else:
        with io.open(filtered_filename, 'w') as out_file:
            for seq in skbio.io.read(fastx, **io_kwargs):
                if exclude_reads:
                    if counter not in filtered_rows:
                        seq.write(out_file, **io_kwargs)
                else:
                    if counter in filtered_rows:
                        seq.write(out_file, **io_kwargs)
                counter += 1
