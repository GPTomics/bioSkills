#!/usr/bin/env python3
'''
cfDNA fragmentomics analysis for cancer detection.
'''

import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def calculate_fragment_metrics(bam_path):
    '''Calculate DELFI-style fragment metrics.'''
    bam = pysam.AlignmentFile(bam_path, 'rb')
    sizes = []

    for read in bam.fetch():
        if read.is_proper_pair and not read.is_secondary and read.template_length > 0:
            sizes.append(read.template_length)

    bam.close()
    sizes = np.array(sizes)

    short = np.sum((sizes >= 100) & (sizes <= 150))
    long = np.sum((sizes >= 151) & (sizes <= 220))

    return {
        'total_fragments': len(sizes),
        'median_size': np.median(sizes),
        'mean_size': np.mean(sizes),
        'short_fragments': short,
        'long_fragments': long,
        'short_long_ratio': short / long if long > 0 else np.nan,
        'mono_peak_fraction': np.sum((sizes >= 150) & (sizes <= 180)) / len(sizes)
    }


def calculate_binned_profile(bam_path, bin_size=5_000_000, chromosomes=None):
    '''Calculate fragment profiles in genomic bins (DELFI-style).'''
    if chromosomes is None:
        chromosomes = [f'chr{i}' for i in range(1, 23)]

    bam = pysam.AlignmentFile(bam_path, 'rb')
    profiles = []

    for chrom in chromosomes:
        try:
            chrom_len = bam.get_reference_length(chrom)
        except Exception:
            continue

        n_bins = (chrom_len // bin_size) + 1

        for bin_idx in range(n_bins):
            start = bin_idx * bin_size
            end = min((bin_idx + 1) * bin_size, chrom_len)

            short_count, long_count = 0, 0

            for read in bam.fetch(chrom, start, end):
                if not read.is_proper_pair or read.is_secondary:
                    continue
                size = read.template_length
                if size <= 0:
                    continue
                if 100 <= size <= 150:
                    short_count += 1
                elif 151 <= size <= 220:
                    long_count += 1

            ratio = short_count / long_count if long_count > 0 else np.nan

            profiles.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'short': short_count,
                'long': long_count,
                'ratio': ratio
            })

    bam.close()
    return pd.DataFrame(profiles)


def plot_genome_profile(profile_df, output_file):
    '''Plot genome-wide fragmentation profile.'''
    fig, ax = plt.subplots(figsize=(15, 5))

    chroms = profile_df['chrom'].unique()
    x_offset = 0
    x_ticks, x_labels = [], []

    for chrom in chroms:
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
    ax.set_ylabel('Short/Long Fragment Ratio')
    ax.set_title('Genome-wide Fragmentation Profile')

    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()


def compare_to_reference(sample_profile, reference_profile):
    '''Compare sample to healthy reference.'''
    merged = sample_profile.merge(
        reference_profile,
        on=['chrom', 'start', 'end'],
        suffixes=('_sample', '_ref')
    )

    merged['z_score'] = (merged['ratio_sample'] - merged['ratio_ref'].mean()) / merged['ratio_ref'].std()
    merged['deviation'] = merged['ratio_sample'] - merged['ratio_ref']

    return merged


if __name__ == '__main__':
    print('Fragment Analysis for Cancer Detection')
    print('=' * 40)
    print('1. calculate_fragment_metrics() - Basic fragment statistics')
    print('2. calculate_binned_profile() - Genome-wide fragmentation')
    print('3. plot_genome_profile() - Visualize profile')
    print('4. compare_to_reference() - Compare to healthy')
