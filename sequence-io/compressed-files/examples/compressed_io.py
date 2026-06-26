#!/usr/bin/env python3
'''Read, write, and index compressed sequence files (gzip / BGZF).'''
# Reference: biopython 1.83+ | Verify API if version differs

import gzip
import os
from Bio import SeqIO, bgzf
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def build_demo_records(n=5):
    return [SeqRecord(Seq('ACGTACGTACGT' * (i + 1)), id=f'gene_{i:03d}', description=f'gene_{i:03d} demo') for i in range(n)]


def write_gzip(records, filepath, fmt='fasta'):
    '''Write records to a plain gzip stream. Text mode is required for SeqIO.'''
    with gzip.open(filepath, 'wt') as handle:
        return SeqIO.write(records, handle, fmt)


def write_bgzf(records, filepath, fmt='fasta'):
    '''Write records to BGZF, the only seekable/indexable compressed format.'''
    with bgzf.open(filepath, 'wt') as handle:
        return SeqIO.write(records, handle, fmt)


def count_compressed(filepath, opener, fmt='fasta'):
    '''Stream-count records without holding them in memory.'''
    with opener(filepath, 'rt') as handle:
        return sum(1 for _ in SeqIO.parse(handle, fmt))


def fetch_from_bgzf(filepath, record_id, fmt='fasta'):
    '''Random access by id via the BGZF offset index.'''
    index = SeqIO.index(filepath, fmt)
    seq = str(index[record_id].seq)
    index.close()
    return seq


def convert_gzip_to_bgzf(gz_path, bgz_path, fmt='fasta'):
    '''Re-compress a plain gzip file as BGZF so it becomes indexable.'''
    with gzip.open(gz_path, 'rt') as in_handle:
        with bgzf.open(bgz_path, 'wt') as out_handle:
            return SeqIO.write(SeqIO.parse(in_handle, fmt), out_handle, fmt)


if __name__ == '__main__':
    records = build_demo_records()
    write_gzip(records, 'demo.fasta.gz')
    write_bgzf(records, 'demo.fasta.bgz')

    print('plain gzip records:', count_compressed('demo.fasta.gz', gzip.open))
    print('BGZF records:', count_compressed('demo.fasta.bgz', bgzf.open))
    print('random access gene_003 from BGZF:', fetch_from_bgzf('demo.fasta.bgz', 'gene_003'))

    voffset = bgzf.make_virtual_offset(100, 7)
    print('virtual offset 100<<16|7 =', voffset, '->', bgzf.split_virtual_offset(voffset))

    # BGZF is seekable by virtual offset; plain gzip is not indexable at all.
    with bgzf.open('demo.fasta.bgz', 'rt') as handle:
        handle.readline()
        saved = handle.tell()
        print('saved virtual offset, re-seek line:', handle.seek(saved) == saved)

    try:
        SeqIO.index('demo.fasta.gz', 'fasta')
    except ValueError as err:
        print('plain gzip index rejected:', str(err))

    convert_gzip_to_bgzf('demo.fasta.gz', 'converted.fasta.bgz')
    print('converted gzip->BGZF, gene_004:', fetch_from_bgzf('converted.fasta.bgz', 'gene_004')[:12], '...')

    for f in ('demo.fasta.gz', 'demo.fasta.bgz', 'converted.fasta.bgz'):
        os.remove(f)
