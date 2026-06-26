#!/usr/bin/env python3
'''Calculate assembly and sequence statistics (N50/L50, auN, GC content).'''
# Reference: biopython 1.83+ | Verify API if version differs

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import gzip
import statistics
import sys
import os
import tempfile

def n50_l50(sorted_lengths):
    '''N50 = minimal length covering >=50% of total; L50 = contig count to reach it.'''
    half = sum(sorted_lengths) / 2
    cumsum = 0
    for count, length in enumerate(sorted_lengths, start=1):
        cumsum += length
        if cumsum >= half:
            return length, count
    return 0, 0

def calculate_stats(filepath, format_type='fasta'):
    '''Compute contiguity and composition stats in a single pass.'''
    opener = gzip.open if str(filepath).endswith('.gz') else open

    lengths = []
    gc_values = []
    n_counts = []
    with opener(filepath, 'rt') as handle:
        for record in SeqIO.parse(handle, format_type):
            lengths.append(len(record.seq))
            gc_values.append(gc_fraction(record.seq, ambiguous='remove'))
            n_counts.append(str(record.seq).upper().count('N'))

    total_bases = sum(lengths)
    sorted_lengths = sorted(lengths, reverse=True)
    n50, l50 = n50_l50(sorted_lengths)
    aun = sum(length * length for length in lengths) / total_bases if total_bases else 0

    return {
        'num_sequences': len(lengths),
        'total_bases': total_bases,
        'min_length': min(lengths),
        'max_length': max(lengths),
        'mean_length': statistics.mean(lengths),
        'median_length': statistics.median(lengths),
        'n50': n50,
        'l50': l50,
        'aun': aun,
        'mean_gc_pct': statistics.mean(gc_values) * 100,
        'total_ns': sum(n_counts),
        'n_fraction_pct': sum(n_counts) / total_bases * 100 if total_bases else 0,
    }

def print_stats(stats):
    print(f'{"Metric":<20} {"Value":>15}')
    print('-' * 36)
    for key, value in stats.items():
        if isinstance(value, float):
            print(f'{key:<20} {value:>15.2f}')
        else:
            print(f'{key:<20} {value:>15,}')

DEMO_CONTIGS = {'contig1': 1500, 'contig2': 800, 'contig3': 800, 'contig4': 400, 'contig5': 200}

if __name__ == '__main__':
    if len(sys.argv) > 1:
        print_stats(calculate_stats(sys.argv[1]))
    else:
        workdir = tempfile.mkdtemp(prefix='seq_stats_')
        demo = os.path.join(workdir, 'sequences.fasta')
        with open(demo, 'w') as fh:
            for name, length in DEMO_CONTIGS.items():
                fh.write(f'>{name}\n' + 'ATGC' * (length // 4) + '\n')
        print_stats(calculate_stats(demo))
        os.remove(demo)
        os.rmdir(workdir)
