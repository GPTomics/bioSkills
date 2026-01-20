#!/usr/bin/env python3
'''Visualize CNV profiles'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

def plot_cnv_profile(cnr_file, output_file='cnv_profile.png'):
    '''Plot genome-wide CNV profile'''
    df = pd.read_csv(cnr_file, sep='\t')

    chrom_order = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    df['chrom'] = pd.Categorical(df['chromosome'], categories=chrom_order, ordered=True)
    df = df.sort_values(['chrom', 'start'])

    # Calculate cumulative position
    chrom_sizes = df.groupby('chromosome')['end'].max()
    chrom_offsets = {}
    offset = 0
    for chrom in chrom_order:
        if chrom in chrom_sizes.index:
            chrom_offsets[chrom] = offset
            offset += chrom_sizes[chrom]

    df['plot_pos'] = df.apply(lambda x: chrom_offsets.get(x['chromosome'], 0) + x['start'], axis=1)

    # Plot
    fig, ax = plt.subplots(figsize=(16, 4))

    colors = ['steelblue' if i % 2 == 0 else 'coral' for i in range(len(chrom_order))]
    for i, chrom in enumerate(chrom_order):
        mask = df['chromosome'] == chrom
        if mask.any():
            ax.scatter(df.loc[mask, 'plot_pos'], df.loc[mask, 'log2'],
                      c=colors[i], s=1, alpha=0.5)

    ax.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax.axhline(0.5, color='red', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axhline(-0.5, color='blue', linestyle='--', linewidth=0.5, alpha=0.5)

    ax.set_xlabel('Genomic Position')
    ax.set_ylabel('Log2 Ratio')
    ax.set_ylim(-2, 2)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f'Saved: {output_file}')

if __name__ == '__main__':
    cnr_file = sys.argv[1] if len(sys.argv) > 1 else 'sample.cnr'
    plot_cnv_profile(cnr_file)
