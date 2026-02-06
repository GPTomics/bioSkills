#!/usr/bin/env python3
'''
Compute repeat statistics from RepeatMasker output.
'''

import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO


def parse_repeatmasker_out(out_file):
    '''Parse RepeatMasker .out file into a DataFrame.

    The .out file has a 3-line header followed by space-delimited columns:
    score, div%, del%, ins%, query, begin, end, left, strand, repeat, class/family, ...
    '''
    records = []
    with open(out_file) as f:
        for i, line in enumerate(f):
            if i < 3:
                continue
            parts = line.split()
            if len(parts) < 15:
                continue
            records.append({
                'score': int(parts[0]),
                'perc_div': float(parts[1]),
                'perc_del': float(parts[2]),
                'perc_ins': float(parts[3]),
                'seqid': parts[4],
                'start': int(parts[5]),
                'end': int(parts[6]),
                'strand': '+' if parts[8] == '+' else '-',
                'repeat_name': parts[9],
                'repeat_class': parts[10].split('/')[0],
                'repeat_family': parts[10],
                'length': int(parts[6]) - int(parts[5]) + 1,
            })
    return pd.DataFrame(records)


def genome_size_from_fasta(fasta_file):
    return sum(len(rec.seq) for rec in SeqIO.parse(fasta_file, 'fasta'))


def repeat_summary(rm_df, genome_size):
    '''Summarize repeat content by class and family.'''
    total_masked = rm_df['length'].sum()

    print(f'=== Repeat Summary ===')
    print(f'Genome size: {genome_size:,} bp')
    print(f'Total masked: {total_masked:,} bp ({total_masked/genome_size:.1%})')
    print(f'Total elements: {len(rm_df):,}')

    class_stats = rm_df.groupby('repeat_class').agg(
        count=('repeat_name', 'count'),
        total_bp=('length', 'sum'),
        mean_div=('perc_div', 'mean'),
    ).sort_values('total_bp', ascending=False)
    class_stats['pct_genome'] = class_stats['total_bp'] / genome_size * 100

    print(f'\nBy class:')
    for cls, row in class_stats.iterrows():
        print(f'  {cls}: {row["count"]:,.0f} elements, {row["total_bp"]:,.0f} bp ({row["pct_genome"]:.1f}%), mean div {row["mean_div"]:.1f}%')

    return class_stats


def repeat_landscape(rm_df, output_file='repeat_landscape.png'):
    '''Plot repeat divergence landscape showing TE accumulation history.

    Low divergence = recent insertions.
    High divergence = ancient insertions.
    Peaks indicate bursts of TE activity.
    '''
    fig, ax = plt.subplots(figsize=(12, 6))

    major_classes = {
        'LINE': '#1f77b4', 'SINE': '#ff7f0e', 'LTR': '#2ca02c',
        'DNA': '#d62728', 'Simple_repeat': '#9467bd',
    }

    for cls, color in major_classes.items():
        subset = rm_df[rm_df['repeat_class'] == cls]
        if len(subset) > 0:
            # Weight by element length to show bp contribution
            ax.hist(subset['perc_div'], bins=50, range=(0, 50), weights=subset['length'],
                    alpha=0.6, label=f'{cls} ({len(subset):,})', color=color)

    ax.set_xlabel('Kimura Divergence from Consensus (%)')
    ax.set_ylabel('Base Pairs')
    ax.set_title('Repeat Landscape')
    ax.legend()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    return fig


def repeat_size_distribution(rm_df, output_file='repeat_sizes.png'):
    '''Plot repeat element size distribution by class.'''
    fig, ax = plt.subplots(figsize=(10, 6))

    major_classes = ['LINE', 'SINE', 'LTR', 'DNA']
    data = [rm_df[rm_df['repeat_class'] == cls]['length'].values for cls in major_classes if cls in rm_df['repeat_class'].values]
    labels = [cls for cls in major_classes if cls in rm_df['repeat_class'].values]

    if data:
        ax.boxplot(data, labels=labels, showfliers=False)
        ax.set_ylabel('Element Length (bp)')
        ax.set_title('Repeat Size Distribution by Class')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print('Usage: repeat_stats.py <repeatmasker.out> <genome.fasta>')
        sys.exit(1)

    rm_df = parse_repeatmasker_out(sys.argv[1])
    genome_size = genome_size_from_fasta(sys.argv[2])

    class_stats = repeat_summary(rm_df, genome_size)
    class_stats.to_csv('repeat_class_summary.tsv', sep='\t')

    repeat_landscape(rm_df)
    repeat_size_distribution(rm_df)

    print('\nOutputs: repeat_class_summary.tsv, repeat_landscape.png, repeat_sizes.png')
