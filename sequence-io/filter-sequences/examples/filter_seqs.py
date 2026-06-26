#!/usr/bin/env python3
'''Stream-filter sequences by length and GC content'''
# Reference: biopython 1.83+ | Verify API if version differs

import sys
import os
import gzip
import tempfile
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def open_maybe_gzip(path, mode):
    '''Open plain or gzipped files in text mode for SeqIO'''
    opener = gzip.open if path.endswith('.gz') else open
    return opener(path, mode)

def passes(record, min_length, min_gc, max_gc, gc_mode):
    '''Length and GC predicate; gc_fraction returns a fraction 0-1'''
    if len(record.seq) < min_length:
        return False
    gc = gc_fraction(record.seq, ambiguous=gc_mode)
    return min_gc <= gc <= max_gc

def filter_stream(input_path, output_path, fmt='fasta', min_length=100, min_gc=0.3, max_gc=0.7, gc_mode='ignore'):
    '''Stream records through the predicate so the whole file never loads'''
    with open_maybe_gzip(input_path, 'rt') as fin, open_maybe_gzip(output_path, 'wt') as fout:
        records = SeqIO.parse(fin, fmt)
        survivors = (rec for rec in records if passes(rec, min_length, min_gc, max_gc, gc_mode))
        count = SeqIO.write(survivors, fout, fmt)
    print(f'Wrote {count} sequences to {output_path}')
    return count

DEMO_FASTA = '''>keep_midGC_long
ATGCGATCGATCGATCGGCGCATATGCGCGCATGCATGCATCGGCATCGATCGATCGATCGATCGATCGATCGATCGCATGCATCGATCGATCGATCGATCG
>drop_short
ATGCGATCGATCG
>keep_softmasked
atgcgatcgatcgatcggcgcatatgcgcgcatgcatgcatcggcatcgatcgatcgatcgatcgatcgatcgatcgcatgcatcgatcgatcgatcgatcg
>drop_lowGC
ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT
'''

if __name__ == '__main__':
    if len(sys.argv) > 1:
        filter_stream(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else 'filtered.fasta')
    else:
        workdir = tempfile.mkdtemp(prefix='filter_seqs_')
        demo_in = os.path.join(workdir, 'input.fasta')
        demo_out = os.path.join(workdir, 'filtered.fasta')
        with open(demo_in, 'w') as fh:
            fh.write(DEMO_FASTA)
        filter_stream(demo_in, demo_out)
        os.remove(demo_in)
        os.remove(demo_out)
        os.rmdir(workdir)
