#!/usr/bin/env python3
'''Calculate sequence statistics'''

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import gzip
import numpy as np

def calculate_stats(filepath, format_type='fasta'):
    '''Calculate comprehensive statistics'''
    opener = gzip.open if str(filepath).endswith('.gz') else open

    lengths = []
    gc_contents = []
    n_counts = []

    with opener(filepath, 'rt') as handle:
        for record in SeqIO.parse(handle, format_type):
            seq = str(record.seq).upper()
            lengths.append(len(seq))
            gc_contents.append(gc_fraction(record.seq) * 100)
            n_counts.append(seq.count('N'))

    lengths = np.array(lengths)
    total_bases = sum(lengths)

    # N50 calculation
    sorted_lengths = np.sort(lengths)[::-1]
    cumsum = np.cumsum(sorted_lengths)
    n50_idx = np.searchsorted(cumsum, total_bases / 2)
    n50 = sorted_lengths[n50_idx]

    return {
        'num_sequences': len(lengths),
        'total_bases': total_bases,
        'min_length': int(np.min(lengths)),
        'max_length': int(np.max(lengths)),
        'mean_length': np.mean(lengths),
        'median_length': np.median(lengths),
        'n50': int(n50),
        'mean_gc': np.mean(gc_contents),
        'total_ns': sum(n_counts),
        'n_fraction': sum(n_counts) / total_bases * 100
    }

def print_stats(stats):
    '''Pretty print statistics'''
    print(f"{'Metric':<20} {'Value':>15}")
    print('-' * 36)
    for key, value in stats.items():
        if isinstance(value, float):
            print(f'{key:<20} {value:>15.2f}')
        else:
            print(f'{key:<20} {value:>15,}')

if __name__ == '__main__':
    import sys
    filepath = sys.argv[1] if len(sys.argv) > 1 else 'sequences.fasta'
    stats = calculate_stats(filepath)
    print_stats(stats)
