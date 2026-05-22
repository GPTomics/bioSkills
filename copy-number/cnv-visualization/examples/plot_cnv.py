#!/usr/bin/env python3
'''Genome-wide CNV profile with a paired B-allele-frequency panel.

A log2-only plot cannot show loss of heterozygosity and misleads if the diploid
baseline is wrong. This script stacks log2 copy ratio over a BAF track so that
copy-neutral LOH (log2 ~ 0 but split BAF) is visible. Pass the caller .cnr/.cns
plus a het-SNP BAF table.
'''
# Reference: matplotlib 3.8+, numpy 1.26+, pandas 2.2+ | Verify API if version differs

import sys
import pandas as pd
import matplotlib.pyplot as plt

CHROMS = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']


def cumulative_offsets(positions_by_chrom):
    '''Build cumulative genomic x-offsets in fixed chromosome order.'''
    offsets, cum = {}, 0
    for c in CHROMS:
        if c in positions_by_chrom and positions_by_chrom[c] > 0:
            offsets[c] = cum
            cum += positions_by_chrom[c]
    return offsets, cum


def plot_cnv_with_baf(cnr_file, cns_file, baf_file=None, output='cnv_profile.png',
                      ploidy_baseline=0.0, gain=0.3, loss=-0.3):
    '''Genome-wide log2 profile with segments, optional BAF panel.

    ploidy_baseline shifts the log2 zero line; set it from the caller ploidy
    estimate rather than the data mode to avoid centering artifacts.
    '''
    cnr = pd.read_csv(cnr_file, sep='\t')
    cns = pd.read_csv(cns_file, sep='\t')
    sizes = cnr.groupby('chromosome')['end'].max().to_dict()
    offsets, genome_len = cumulative_offsets(sizes)

    cnr = cnr[cnr['chromosome'].isin(offsets)].copy()
    cnr['x'] = cnr.apply(lambda r: offsets[r['chromosome']] + r['start'], axis=1)

    n_panels = 2 if baf_file else 1
    fig, axes = plt.subplots(n_panels, 1, figsize=(16, 3 * n_panels + 1), sharex=True,
                             squeeze=False)
    ax_cn = axes[0][0]

    ax_cn.scatter(cnr['x'], cnr['log2'] - ploidy_baseline, s=1, c='0.7', alpha=0.5,
                  rasterized=True)
    for _, seg in cns.iterrows():
        if seg['chromosome'] not in offsets:
            continue
        x0 = offsets[seg['chromosome']] + seg['start']
        x1 = offsets[seg['chromosome']] + seg['end']
        lr = seg['log2'] - ploidy_baseline
        color = 'red' if lr > gain else 'blue' if lr < loss else '0.3'
        ax_cn.hlines(lr, x0, x1, colors=color, linewidth=2.5)
    ax_cn.axhline(0, color='black', linewidth=0.6)
    ax_cn.set_ylim(-2, 2)
    ax_cn.set_ylabel('log2 copy ratio')

    if baf_file:
        baf = pd.read_csv(baf_file, sep='\t')
        baf = baf[baf['chromosome'].isin(offsets)].copy()
        baf['x'] = baf.apply(lambda r: offsets[r['chromosome']] + r['position'], axis=1)
        ax_baf = axes[1][0]
        # Plot BAF and its mirror; split bands away from 0.5 indicate LOH / imbalance.
        ax_baf.scatter(baf['x'], baf['baf'], s=2, c='0.4', alpha=0.5, rasterized=True)
        ax_baf.scatter(baf['x'], 1 - baf['baf'], s=2, c='0.4', alpha=0.5, rasterized=True)
        ax_baf.axhline(0.5, color='black', linewidth=0.6)
        ax_baf.set_ylim(0, 1)
        ax_baf.set_ylabel('B-allele frequency')

    for c, x in offsets.items():
        for row in axes:
            row[0].axvline(x, color='0.9', linewidth=0.5)
    axes[-1][0].set_xlabel('genomic position')
    axes[-1][0].set_xlim(0, genome_len)
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    print(f'Saved: {output}')
    return fig


if __name__ == '__main__':
    cnr = sys.argv[1] if len(sys.argv) > 1 else 'sample.cnr'
    cns = sys.argv[2] if len(sys.argv) > 2 else 'sample.cns'
    baf = sys.argv[3] if len(sys.argv) > 3 else None
    plot_cnv_with_baf(cnr, cns, baf)
