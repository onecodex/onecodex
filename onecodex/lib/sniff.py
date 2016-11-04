#!/usr/bin/env python
"""
Sniffing functions for FASTA/FASTQ files
"""
from __future__ import print_function, division

import gzip
import os
import re
from collections import Counter

COMMON_NA = set('ACGNTUX')
IUPAC_NA = set('ABCDGHIKMNRSTUVWXY')
IUPAC_AA = set('ABCDEFGHIKLMNPQRSTUVWXYZ*')


def sniff_file(filename, compress=None):
    """
    Given a sequencing file, return a JSON blob of relevant summary statistics
    about that file.
    """
    if not os.path.exists(filename):
        return {'file_type': 'bad', 'msg': 'File does not exist'}
    elif os.path.getsize(filename) < 35:
        return {'file_type': 'bad', 'msg': 'File is too small'}

    if compress is None:
        with open(filename, 'r') as seq_file:
            start = seq_file.read(1)
            data = seq_file.read(1000000)
    elif compress == 'gzip':
        with gzip.open(filename, 'r') as seq_file:
            start = seq_file.read(1)
            data = seq_file.read(1000000)

    # TODO: check for bzip?
    if start == '\x1f' and compress is None:
        # it was a gzip file, try opening it that way
        return sniff_file(filename, 'gzip')

    status = sniff(start, data)

    if compress == 'gzip':
        status['compression'] = 'gzip'
    else:
        status['compression'] = 'none'
    return status


def sniff(start, data):
    """
    Given the first byte (start) and an unspecified (maybe not all) amount of
    a FASTA or FASTQ, return summary statistics.
    """
    # scan through the file and get ids/seq_counts (and quality info)
    if start == '>':
        seq_count, ids, status = read_fasta(data)
    elif start == '@':
        seq_count, ids, status = read_fastq(data)
    else:
        return {'file_type': 'bad', 'msg': 'File is not a valid FASTA or FASTQ file'}

    num_recs = len(ids)
    if num_recs < 1:
        return {'file_type': 'bad', 'msg': 'No records found in file'}
    elif sum(seq_count.values()) < 1:
        return {'file_type': 'bad', 'msg': 'No sequence data found in file'}
    status.update(sniff_bases(seq_count, num_recs))
    status.update(sniff_ids(ids))

    return status


def sniff_bases(seq_count, num_recs):
    """
    Examine the basepair Counter statistics from a reads file and return
    pertinent statistics.
    """
    status = {}

    # check if there are newlines in the sequence
    status['seq_multiline'] = seq_count.pop('\n', 0) > 0 or seq_count.pop('\r', 0) > 0

    # TODO: strip other whitespace?
    status['seq_est_avg_len'] = sum(seq_count.values()) / num_recs

    # strip out gaps (and record if they were there)
    status['seq_has_gaps'] = seq_count.pop('.', 0) + seq_count.pop('-', 0) > 0

    # check if lowercase letters are present in the sequence; uppercase them if so
    status['seq_has_lowercase'] = False
    for k in seq_count.keys():
        if k.islower():
            status['seq_has_lowercase'] = True
            seq_count[k.upper()] += seq_count.pop(k)

    seq_set = set(seq_count)

    # determine the type of the sequence based on the most abundant bases
    seq_common = set(i[0] for i in seq_count.most_common(5))
    if seq_common.issubset(COMMON_NA):
        status['seq_type'] = 'dna' if 'U' not in seq_common else 'rna'
    else:
        status['seq_type'] = 'aa'

    if status['seq_type'] in ['dna', 'rna']:
        # do a gc % calculation (accounting for IUPAC codes)
        gc = seq_count['G'] + seq_count['C'] + seq_count['S']
        gc += 0.5 * (seq_count['K'] + seq_count['M'] + seq_count['R'] + seq_count['Y'] +
                     seq_count['N'] + seq_count['X'])
        gc += 0.333 * (seq_count['D'] + seq_count['H']) + 0.667 * (seq_count['B'] + seq_count['V'])

        # check for IUPAC nucleotides
        if IUPAC_NA.difference(COMMON_NA).isdisjoint(seq_set):
            # doesn't have any of the "special" IUPAC letters
            status['seq_has_iupac'] = False
            status['seq_est_gc'] = gc / sum(seq_count.values())
        elif seq_set.issubset(IUPAC_NA):
            # is a strict subset of the IUPAC letters
            status['seq_has_iupac'] = True
            status['seq_est_gc'] = gc / sum(seq_count.values())
        else:
            # if it has non-nucleotide codes in it, it's not dna/rna
            status['seq_type'] = 'aa'

    # make sure that AA files doesn't have weird non-IUPAC characters
    if status['seq_type'] == 'aa':
        # there's a J, O, or some random character in the sequences
        status['seq_has_noniupac'] = not seq_set.issubset(IUPAC_AA)

    status['seq_has_unknowns'] = ('N' in seq_set and status['seq_type'] != 'aa') or 'X' in seq_set

    return status


def sniff_ids(ids):
    """
    Scan through the IDs from a FASTA/Q and return pertinant statistics.

    """
    status = {}

    # check for interleaving (replace 2 with 1 and see if every two are duplicates)
    # this won't catch a single unpaired read at the end of the file (because we don't know if the
    # file is longer than 1 Mb and we don't want to miscall because one half of a read was cut out)
    singled = [i.replace('2', '1') for i in ids]
    status['interleaved'] = all(singled[2 * i] == singled[2 * i + 1] for
                                i in range(len(singled) // 2)) and len(ids) > 1
    return status

    # TODO: return id_est_len to help estimate # of sequences in file
    # TODO: determine id type? (sequencer, assembler, database ...)

    # # https://en.wikipedia.org/wiki/FASTA_format#Sequence_identifiers
    # db_starters = {  # all ids start with one of the following like `id|...`
    #     'gb': 'GenBank',
    #     'emb': 'EMBL Data Library',
    #     'dbj': 'DNA Database of Japan',
    #     'pir': 'NBRF PIR',
    #     'sp': 'SWISS-PROT',
    #     'pdb': 'Brookhaven Protein Data Bank',
    #     'pat': 'Patents',
    #     'bbs': 'GenInfo Backbone',
    #     'ref': 'NCBI Reference',
    #     'lcl': 'Local Sequence',  # ?
    #     'tr': 'TrEMBL',
    #     'gnl': 'General database identifier',  # ?
    #     'gi': 'NCBI'  # ?
    # }

    #    assemblers = {
    #        'trinity': r'c\d+_g\d+_\w+',  # http://seqanswers.com/forums/showthread.php?t=43749
    #        'idba': '',
    #        # based off sequences from http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/
    #        'abyss_contig': r'[\d.]+ \d+ \d+ (?P<seq_len>\d+) [f]',
    #        'abyss_scaffold': r'\d+ \d+ (?P<seq_len>\d+)',
    #        'allpaths_contig': r'contig_\d+',
    #        'allpaths_scaffold': r'scaffold_\d+',
    #        'bambus_contig': 'scf\d+_\d+',
    #        'cabog,msrca_contig': 'ctg\d+',
    #        'bambus,cabog,msrca_scaffold': 'scf\d+',
    #        'sga_contig': 'contig-\d+ (?P<seq_len>\d+) \d',
    #        'sga_scaffold': 'scaffold-\d+ (?P<seq_len>\d+)',
    #        'soapdenovo_contig': 'C[\d.]+ C\d+ \d (?P<seq_len>\d+) [f]',
    #        'soapdenovo_scaffold': 'C\d+ [\d.]+',
    #        'velvet_contig': 'velvet.\d+.\d+ \d+ (?P<seq_len>\d+)',
    #        'velvet_scaffold': 'velvet.\d+ NODE_\d+_length_(?P<seq_len>\d+)_cov_(?P<seq_cov>[\d.]+)'  # noqa
    #        # based off sequences from http://gigadb.org/dataset/100060
    #
    #        # see badger for 454 sequencing files
    #
    #        # based off http://www.ncbi.nlm.nih.gov/sra?term=SRA026860
    #        # http://korflab.ucdavis.edu/Datasets/Assemblathon/Assemblathon1/Entries/
    #    }

    # sequencers = {
    #     # http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/
    #     #        Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm
    #     # http://wiki.christophchamp.com/index.php/FASTQ_format
    #     'illumina': '',
    # }


def read_fasta(data):
    fasta_re = re.compile(r"""
        (?P<id>[^\n]+)\n  # the identifier line
        (?P<seq>[^>]+)  # the sequence
        (?:\n>|\Z)  # start of the next record
    """, re.VERBOSE)
    ids = []
    seq_count = Counter()
    for match in fasta_re.finditer(data):
        rec = match.groupdict()
        ids.append(rec['id'])
        seq_count.update(Counter(rec['seq'].rstrip()))
    return seq_count, ids, {'file_type': 'fasta'}


def read_fastq(data):
    fastq_re = re.compile(r"""
        (?P<id>[^\n]+)\n
        (?P<seq>[^\n]+)\n
        \+(?P<id2>[^\n]*)\n
        (?P<qual>[^\n]+)
        (?:\n@|\Z)
    """, re.DOTALL + re.VERBOSE)
    status = {'file_type': 'fastq'}
    qual_set = set()
    ids = []
    seq_count = Counter()
    for match in fastq_re.finditer(data):
        rec = match.groupdict()
        if rec['id'] != rec['id2']:
            # once the qual_ids don't match, we always report `nonmatch`
            if rec['id2'] == '' and status.get('qual_ids') != 'nonmatch':
                status['qual_ids'] = 'blank_second'
            else:
                status['qual_ids'] = 'nonmatch'
        ids.append(rec['id'])
        qual_set.update(rec['qual'])
        seq_count.update(Counter(rec['seq']))

    if 'qual_ids' not in status:
        status['qual_ids'] = 'match'

    qual_set = qual_set.difference({'\n', '\r', ' ', '\t'})
    # https://en.wikipedia.org/wiki/FASTQ_format#Encoding
    printable = [chr(i) for i in range(33, 127)]
    if qual_set.issubset(printable[0:41]):
        # we call sanger before illumina 1.8 b/c it's a technically more restricted subset
        status['qual_type'] = 'sanger'
    elif qual_set.issubset(printable[0:42]):
        status['qual_type'] = 'illumina 1.8'
    elif qual_set.issubset(printable[33:72]):
        status['qual_type'] = 'illumina 1.5'
    elif qual_set.issubset(printable[31:72]):
        status['qual_type'] = 'illumina 1.3'
    elif qual_set.issubset(printable[26:72]):
        status['qual_type'] = 'solexa'
    else:
        status['qual_type'] = 'bad'

    return seq_count, ids, status


if __name__ == '__main__':
    import argparse
    import json

    parser = argparse.ArgumentParser(description='Sniff sequencing files for information.')
    parser.add_argument('filename', default=None, nargs='?', help='File to sniff')

    args = parser.parse_args()
    print(json.dumps(sniff_file(args.filename)))
