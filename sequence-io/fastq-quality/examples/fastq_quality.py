#!/usr/bin/env python3
'''Analyze FASTQ quality scores'''

from Bio import SeqIO
import numpy as np
import gzip

def parse_fastq(filepath):
    '''Parse FASTQ, handling gzip'''
    opener = gzip.open if str(filepath).endswith('.gz') else open
    with opener(filepath, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fastq'):
            yield record

def calculate_quality_stats(filepath, sample_size=10000):
    '''Calculate quality statistics'''
    qualities = []
    lengths = []
    gc_contents = []

    for i, record in enumerate(parse_fastq(filepath)):
        if i >= sample_size:
            break
        quals = record.letter_annotations['phred_quality']
        qualities.append(np.mean(quals))
        lengths.append(len(record))
        gc = (record.seq.count('G') + record.seq.count('C')) / len(record) * 100
        gc_contents.append(gc)

    return {
        'n_reads': i + 1,
        'mean_quality': np.mean(qualities),
        'mean_length': np.mean(lengths),
        'mean_gc': np.mean(gc_contents),
        'q30_fraction': np.mean([q >= 30 for q in qualities])
    }

def per_position_quality(filepath, sample_size=10000):
    '''Calculate per-position quality'''
    max_len = 0
    position_quals = {}

    for i, record in enumerate(parse_fastq(filepath)):
        if i >= sample_size:
            break
        quals = record.letter_annotations['phred_quality']
        for pos, q in enumerate(quals):
            if pos not in position_quals:
                position_quals[pos] = []
            position_quals[pos].append(q)

    positions = sorted(position_quals.keys())
    means = [np.mean(position_quals[p]) for p in positions]
    return positions, means

if __name__ == '__main__':
    import sys
    filepath = sys.argv[1] if len(sys.argv) > 1 else 'reads.fastq.gz'
    stats = calculate_quality_stats(filepath)
    print(f'Quality stats for {filepath}:')
    for k, v in stats.items():
        print(f'  {k}: {v:.2f}' if isinstance(v, float) else f'  {k}: {v}')
