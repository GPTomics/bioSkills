#!/usr/bin/env python
# Reference: pysam 0.22+, pybedtools 0.9+, numpy 1.26+ | Verify API if version differs
# ChIP-seq QC battery: FRiP, NRF/PBC1/PBC2 (pre-dedup), fragment-size diagnostic,
# hyper-ChIPable artifact flagging. Run on the pre-MarkDuplicates BAM for
# library-complexity metrics; post-filtering BAM for FRiP.

import argparse
import sys
import pysam
import pybedtools
import numpy as np


def calculate_frip(bam_file, peak_file):
    '''FRiP = reads in peaks / total mapped primary reads. ENCODE thresholds:
    > 1% (minimum, TF), > 5% (good TF), > 10-20% (histone, broad), > 15% (H3K4me3).
    Compute on the post-filtering, blacklist-cleaned peak file.'''
    bam = pysam.AlignmentFile(bam_file, 'rb')
    total_reads = bam.count(read_callback=lambda r: not r.is_unmapped and not r.is_secondary)

    peaks = pybedtools.BedTool(peak_file)
    reads_in_peaks = 0
    for peak in peaks:
        reads_in_peaks += bam.count(peak.chrom, peak.start, peak.end)

    bam.close()
    return reads_in_peaks / total_reads if total_reads > 0 else 0


def calculate_complexity(bam_file):
    '''Library complexity: NRF, PBC1, PBC2. Per ENCODE, compute on the
    PRE-MarkDuplicates BAM (post-MarkDuplicates BAMs report NRF ~ 1.0 by
    construction, which is uninformative). NRF > 0.8 ideal, < 0.5 severe;
    PBC1 > 0.8 ideal; PBC2 > 3 ideal.'''
    bam = pysam.AlignmentFile(bam_file, 'rb')
    position_counts = {}
    total = 0
    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        total += 1
        pos = (read.reference_name, read.reference_start, read.is_reverse)
        position_counts[pos] = position_counts.get(pos, 0) + 1
    bam.close()

    m1 = sum(1 for v in position_counts.values() if v == 1)
    m2 = sum(1 for v in position_counts.values() if v == 2)
    m_distinct = len(position_counts)

    nrf = m_distinct / total if total > 0 else 0
    pbc1 = m1 / m_distinct if m_distinct > 0 else 0
    pbc2 = m1 / m2 if m2 > 0 else float('inf')

    return {'NRF': nrf, 'PBC1': pbc1, 'PBC2': pbc2,
            'M1': m1, 'M2': m2, 'Mdistinct': m_distinct, 'total': total}


def fragment_size_diagnostic(bam_file, max_size=1000):
    '''Fragment-size diagnostic from properly-paired reads. Reports the mode,
    fraction sub-nucleosomal (<100 bp; expected for TF ChIP), and fraction
    nucleosomal (140-180 bp; expected for histone ChIP). Flat distribution
    indicates over-sonication; biology cannot be rescued.'''
    bam = pysam.AlignmentFile(bam_file, 'rb')
    sizes = []
    for read in bam.fetch():
        if read.is_proper_pair and read.template_length > 0:
            tlen = read.template_length
            if 0 < tlen <= max_size:
                sizes.append(tlen)
    bam.close()

    sizes_arr = np.array(sizes)
    if len(sizes_arr) == 0:
        return None

    sub_nuc = float(np.mean(sizes_arr < 100))
    mono_nuc = float(np.mean((sizes_arr >= 140) & (sizes_arr <= 180)))
    over_300 = float(np.mean(sizes_arr > 300))
    median_size = float(np.median(sizes_arr))

    hist, edges = np.histogram(sizes_arr, bins=range(0, max_size + 1, 10))
    mode_size = int(edges[np.argmax(hist)])

    classification = 'unknown'
    if sub_nuc > 0.30 and mono_nuc < 0.20:
        classification = 'TF-pattern (sub-nucleosomal dominant)'
    elif mono_nuc > 0.20 and sub_nuc < 0.15:
        classification = 'histone-pattern (nucleosomal dominant)'
    elif sub_nuc < 0.05 and mono_nuc < 0.10 and over_300 > 0.40:
        classification = 'over-sonicated (flat distribution; unrecoverable)'
    elif sub_nuc < 0.10 and mono_nuc > 0.15:
        classification = 'mostly nucleosomal (TF? check for trapping)'

    return {'median_size': median_size, 'mode_size': mode_size,
            'frac_sub_nucleosomal_lt100': sub_nuc,
            'frac_mono_nucleosomal_140_180': mono_nuc,
            'frac_over_300': over_300, 'n_pairs': len(sizes_arr),
            'classification': classification}


def hyper_chipable_check(peak_file, input_bigwig_summary_tsv, top_pct=1.0):
    '''Flag peaks falling into top-N% input signal regions (custom hyper-ChIPable
    blacklist). Teytelman 2013 / Park 2013: rRNA, tRNA, histone clusters, mtDNA,
    abundant housekeeping genes appear "bound" even by untagged GFP. ENCODE
    blacklist v2 misses these; cell-type-specific input-based blacklist catches them.
    `input_bigwig_summary_tsv` is multiBigwigSummary --outRawCounts output.'''
    bg_signal = []
    with open(input_bigwig_summary_tsv) as f:
        header = f.readline()
        for line in f:
            fields = line.rstrip().split('\t')
            try:
                signal = float(fields[3])
                bg_signal.append((fields[0], int(fields[1]), int(fields[2]), signal))
            except (ValueError, IndexError):
                continue

    if not bg_signal:
        return None

    bg_signal.sort(key=lambda x: x[3], reverse=True)
    n_top = max(1, int(len(bg_signal) * top_pct / 100))
    top_regions = bg_signal[:n_top]

    blacklist_bed = '\n'.join(f"{c}\t{s}\t{e}" for c, s, e, _ in top_regions)
    blacklist = pybedtools.BedTool(blacklist_bed, from_string=True)
    peaks = pybedtools.BedTool(peak_file)

    suspicious = peaks.intersect(blacklist, u=True)
    n_total = sum(1 for _ in peaks)
    n_suspicious = sum(1 for _ in suspicious)

    return {'total_peaks': n_total, 'suspicious_peaks': n_suspicious,
            'fraction_suspicious': n_suspicious / n_total if n_total > 0 else 0,
            'top_pct_threshold': top_pct}


def grade(metrics, mark_type='tf'):
    '''Grade metrics against ENCODE thresholds. mark_type: tf | sharp_histone |
    broad_histone. Returns pass/caution/reject per metric.'''
    grades = {}
    if mark_type == 'tf':
        thresholds = {'FRiP': (0.05, 0.01), 'NRF': (0.8, 0.5), 'PBC1': (0.8, 0.5)}
    elif mark_type == 'sharp_histone':
        thresholds = {'FRiP': (0.15, 0.05), 'NRF': (0.8, 0.5), 'PBC1': (0.8, 0.5)}
    elif mark_type == 'broad_histone':
        thresholds = {'FRiP': (0.20, 0.10), 'NRF': (0.8, 0.5), 'PBC1': (0.8, 0.5)}
    else:
        return None

    for metric, (good, marginal) in thresholds.items():
        if metric in metrics:
            v = metrics[metric]
            if v >= good:
                grades[metric] = 'pass'
            elif v >= marginal:
                grades[metric] = 'caution'
            else:
                grades[metric] = 'reject'
    return grades


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ChIP-seq QC battery')
    parser.add_argument('bam', help='Post-filtering BAM for FRiP (deduplicated, blacklist-cleaned)')
    parser.add_argument('--pre-dedup-bam', help='Pre-deduplication BAM for NRF/PBC (required for accurate complexity)')
    parser.add_argument('--peaks', help='Peak file (narrowPeak/broadPeak/bed) for FRiP')
    parser.add_argument('--mark-type', choices=['tf', 'sharp_histone', 'broad_histone'],
                        default='tf', help='Mark category for ENCODE threshold grading')
    parser.add_argument('--hyper-chipable-input', help='multiBigwigSummary --outRawCounts TSV from input.bw')
    args = parser.parse_args()

    metrics = {}
    print(f'BAM: {args.bam}')
    print(f'Mark type: {args.mark_type}')

    fs = fragment_size_diagnostic(args.bam)
    if fs:
        print(f'\nFragment size: median={fs["median_size"]:.0f} bp, '
              f'mode={fs["mode_size"]} bp, n={fs["n_pairs"]}')
        print(f'  sub-nucleosomal (<100): {fs["frac_sub_nucleosomal_lt100"]:.2%}')
        print(f'  mono-nucleosomal (140-180): {fs["frac_mono_nucleosomal_140_180"]:.2%}')
        print(f'  > 300 bp: {fs["frac_over_300"]:.2%}')
        print(f'  classification: {fs["classification"]}')

    if args.peaks:
        frip = calculate_frip(args.bam, args.peaks)
        metrics['FRiP'] = frip
        print(f'\nFRiP: {frip:.4f}')

    if args.pre_dedup_bam:
        complexity = calculate_complexity(args.pre_dedup_bam)
        metrics.update({'NRF': complexity['NRF'], 'PBC1': complexity['PBC1'], 'PBC2': complexity['PBC2']})
        print(f'\nLibrary complexity (pre-dedup):')
        for k, v in complexity.items():
            print(f'  {k}: {v:.4f}' if isinstance(v, float) else f'  {k}: {v}')

    if args.hyper_chipable_input and args.peaks:
        hc = hyper_chipable_check(args.peaks, args.hyper_chipable_input)
        if hc:
            print(f'\nHyper-ChIPable artifact check:')
            print(f'  {hc["suspicious_peaks"]} / {hc["total_peaks"]} peaks in top-{hc["top_pct_threshold"]:.0f}% input signal')
            print(f'  Fraction: {hc["fraction_suspicious"]:.2%}')
            if hc['fraction_suspicious'] > 0.10:
                print(f'  WARNING: high fraction of peaks at hyper-ChIPable regions; consider custom blacklist')

    grades = grade(metrics, args.mark_type)
    if grades:
        print(f'\nENCODE grade ({args.mark_type}):')
        for metric, g in grades.items():
            print(f'  {metric}: {g}')
        if 'reject' in grades.values():
            print('\nOVERALL: REJECT (at least one metric failed)', file=sys.stderr)
            sys.exit(1)
