#!/usr/bin/env python3
'''
cfDNA fragmentomics: custom DELFI-style binned short/long ratio + size summary.

This is a dependency-light, illustrative readout. It is NOT GC-corrected, so the
resulting ratios are not comparable across samples, batches, or library protocols.
For publication-grade features use FinaleToolkit (delfi/end-motifs/wps), whose
delfi command GC-corrects, or the Griffin Snakemake pipeline.
'''
# Reference: numpy 1.26+, pandas 2.2+, pysam 0.22+ | Verify API if version differs

import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SHORT_LO, SHORT_HI = 100, 150          # DELFI short window (Cristiano 2019 Nature 570:385)
LONG_LO, LONG_HI = 151, 220            # DELFI long window
MONO_LO, MONO_HI = 150, 180            # ~167 bp mononucleosome mode (Snyder 2016 Cell 164:57)
DEFAULT_BIN = 5_000_000                # ~5 Mb DELFI bin (Cristiano 2019)


def fragment_metrics(bam_path):
    '''Whole-sample size summary plus the global short/long ratio.'''
    bam = pysam.AlignmentFile(bam_path, 'rb')
    sizes = [r.template_length for r in bam.fetch()
             if r.is_proper_pair and not r.is_secondary and r.template_length > 0]
    bam.close()
    sizes = np.array(sizes)
    short = np.sum((sizes >= SHORT_LO) & (sizes <= SHORT_HI))
    long = np.sum((sizes >= LONG_LO) & (sizes <= LONG_HI))
    return {'total_fragments': len(sizes), 'median_size': np.median(sizes), 'mean_size': np.mean(sizes),
            'short_fragments': int(short), 'long_fragments': int(long),
            'short_long_ratio': short / long if long > 0 else np.nan,
            'mono_peak_fraction': np.sum((sizes >= MONO_LO) & (sizes <= MONO_HI)) / len(sizes)}


def binned_short_long_ratio(bam_path, bin_size=DEFAULT_BIN, chroms=None):
    '''Per-bin short/long ratio across the genome. Not GC-corrected -- illustrative only.'''
    chroms = chroms or [f'chr{i}' for i in range(1, 23)]
    bam = pysam.AlignmentFile(bam_path, 'rb')
    rows = []
    for chrom in chroms:
        if chrom not in bam.references:
            continue
        chrom_len = bam.get_reference_length(chrom)
        n_bins = chrom_len // bin_size + 1
        short = np.zeros(n_bins)
        long = np.zeros(n_bins)
        for read in bam.fetch(chrom):
            if not read.is_proper_pair or read.is_secondary or read.template_length <= 0:
                continue
            size = read.template_length
            b = read.reference_start // bin_size
            if SHORT_LO <= size <= SHORT_HI:
                short[b] += 1
            elif LONG_LO <= size <= LONG_HI:
                long[b] += 1
        ratio = np.divide(short, long, out=np.full(n_bins, np.nan), where=long > 0)
        for i in range(n_bins):
            rows.append({'chrom': chrom, 'start': i * bin_size,
                         'end': min((i + 1) * bin_size, chrom_len),
                         'short': short[i], 'long': long[i], 'ratio': ratio[i]})
    bam.close()
    return pd.DataFrame(rows)


def size_selected_metrics(bam_path, lo=90, hi=150):
    '''In-silico size selection (90-150 bp) to enrich tumor fraction (Mouliere 2018 STM 10:eaat4921).'''
    bam = pysam.AlignmentFile(bam_path, 'rb')
    kept = sum(1 for r in bam.fetch()
               if r.is_proper_pair and not r.is_secondary and lo <= r.template_length <= hi)
    total = sum(1 for r in bam.fetch()
                if r.is_proper_pair and not r.is_secondary and r.template_length > 0)
    bam.close()
    return {'kept_fragments': kept, 'total_fragments': total,
            'retained_fraction': kept / total if total > 0 else np.nan}


def compare_to_reference(sample_profile, reference_profile):
    '''Z-score sample bins against a co-processed healthy reference. Only valid if GC-corrected upstream.'''
    merged = sample_profile.merge(reference_profile, on=['chrom', 'start', 'end'], suffixes=('_sample', '_ref'))
    ref_mean, ref_std = merged['ratio_ref'].mean(), merged['ratio_ref'].std()
    merged['z_score'] = (merged['ratio_sample'] - ref_mean) / ref_std
    merged['deviation'] = merged['ratio_sample'] - merged['ratio_ref']
    return merged


def plot_genome_profile(profile_df, output_file):
    '''Scatter the per-bin short/long ratio along the genome.'''
    fig, ax = plt.subplots(figsize=(15, 5))
    x_offset, x_ticks, x_labels = 0, [], []
    for chrom in profile_df['chrom'].unique():
        chrom_data = profile_df[profile_df['chrom'] == chrom]
        x = np.arange(len(chrom_data)) + x_offset
        ax.scatter(x, chrom_data['ratio'], s=10, alpha=0.7)
        x_ticks.append(x_offset + len(chrom_data) / 2)
        x_labels.append(chrom.replace('chr', ''))
        x_offset += len(chrom_data)
        ax.axvline(x=x_offset, color='gray', alpha=0.3, linewidth=0.5)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Short/Long Fragment Ratio (uncorrected)')
    ax.set_title('Genome-wide Fragmentation Profile')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    print('cfDNA fragmentomics (illustrative, NOT GC-corrected)')
    print('fragment_metrics() / binned_short_long_ratio() / size_selected_metrics() / compare_to_reference()')
    print('For GC-corrected features use FinaleToolkit (delfi/end-motifs/wps) or the Griffin pipeline.')
