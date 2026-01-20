#!/usr/bin/env python3
'''Handle paired-end FASTQ files'''

from Bio import SeqIO
import gzip
from pathlib import Path

def parse_paired_fastq(r1_path, r2_path):
    '''Parse paired FASTQ files together'''
    open_r1 = gzip.open if str(r1_path).endswith('.gz') else open
    open_r2 = gzip.open if str(r2_path).endswith('.gz') else open

    with open_r1(r1_path, 'rt') as f1, open_r2(r2_path, 'rt') as f2:
        for rec1, rec2 in zip(SeqIO.parse(f1, 'fastq'), SeqIO.parse(f2, 'fastq')):
            yield rec1, rec2

def validate_pairs(r1_path, r2_path):
    '''Validate paired-end files have matching IDs'''
    mismatches = 0
    total = 0

    for rec1, rec2 in parse_paired_fastq(r1_path, r2_path):
        total += 1
        id1 = rec1.id.split('/')[0].split(' ')[0]
        id2 = rec2.id.split('/')[0].split(' ')[0]
        if id1 != id2:
            mismatches += 1
            if mismatches <= 5:
                print(f'Mismatch: {rec1.id} vs {rec2.id}')

    print(f'Total pairs: {total}')
    print(f'Mismatches: {mismatches}')
    return mismatches == 0

def interleave_fastq(r1_path, r2_path, output_path):
    '''Interleave paired FASTQ files'''
    opener = gzip.open if str(output_path).endswith('.gz') else open

    count = 0
    with opener(output_path, 'wt') as out:
        for rec1, rec2 in parse_paired_fastq(r1_path, r2_path):
            SeqIO.write([rec1, rec2], out, 'fastq')
            count += 1

    print(f'Wrote {count} pairs to {output_path}')

def deinterleave_fastq(interleaved_path, r1_path, r2_path):
    '''Split interleaved FASTQ to paired files'''
    open_in = gzip.open if str(interleaved_path).endswith('.gz') else open
    open_r1 = gzip.open if str(r1_path).endswith('.gz') else open
    open_r2 = gzip.open if str(r2_path).endswith('.gz') else open

    with open_in(interleaved_path, 'rt') as fin, \
         open_r1(r1_path, 'wt') as f1, open_r2(r2_path, 'wt') as f2:
        records = SeqIO.parse(fin, 'fastq')
        for rec1, rec2 in zip(records, records):
            SeqIO.write(rec1, f1, 'fastq')
            SeqIO.write(rec2, f2, 'fastq')

if __name__ == '__main__':
    import sys
    if len(sys.argv) >= 3:
        validate_pairs(sys.argv[1], sys.argv[2])
