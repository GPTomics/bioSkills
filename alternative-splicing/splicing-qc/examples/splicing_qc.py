#!/usr/bin/env python3
'''
Quality control for splicing analysis.
Assesses junction saturation, coverage, and splice site strength.
'''

import subprocess
import pysam
import pandas as pd
from pathlib import Path
from collections import Counter


def run_junction_saturation(bam_file, bed_file, output_prefix):
    '''
    Run RSeQC junction saturation analysis.
    Checks if sequencing depth is sufficient for junction detection.

    Note: -s flag removed in RSeQC v3.0
    '''
    subprocess.run([
        'junction_saturation.py',
        '-i', bam_file,
        '-r', bed_file,
        '-o', output_prefix
    ], check=True)

    print(f'Junction saturation output: {output_prefix}.junctionSaturation_plot.pdf')
    print('Check if curves plateau - plateau indicates sufficient depth')


def run_junction_annotation(bam_file, bed_file, output_prefix):
    '''
    Classify junctions as known, partial novel, or complete novel.
    High proportion of known junctions indicates good alignment.
    '''
    subprocess.run([
        'junction_annotation.py',
        '-i', bam_file,
        '-r', bed_file,
        '-o', output_prefix
    ], check=True)

    # Parse results
    stats_file = f'{output_prefix}.junction.xls'
    if Path(stats_file).exists():
        stats = pd.read_csv(stats_file, sep='\t')
        return stats


def count_junction_reads(bam_path, min_overhang=8):
    '''
    Count reads supporting each splice junction.

    Args:
        bam_path: Path to BAM file
        min_overhang: Minimum overhang on each side of junction
    '''
    bam = pysam.AlignmentFile(bam_path, 'rb')
    junction_counts = Counter()

    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary:
            continue

        ref_pos = read.reference_start
        query_pos = 0

        for op, length in read.cigartuples:
            if op == 3:  # N = splice junction
                junction = (read.reference_name, ref_pos, ref_pos + length)
                junction_counts[junction] += 1
            if op in [0, 2, 3]:  # M, D, N consume reference
                ref_pos += length
            if op in [0, 1, 4]:  # M, I, S consume query
                query_pos += length

    bam.close()
    return junction_counts


def analyze_junction_coverage(junction_counts):
    '''
    Analyze distribution of junction read counts.
    '''
    counts = list(junction_counts.values())

    stats = {
        'total_junctions': len(counts),
        'junctions_ge_10': sum(1 for c in counts if c >= 10),
        'junctions_ge_20': sum(1 for c in counts if c >= 20),
        'junctions_ge_50': sum(1 for c in counts if c >= 50),
        'median_reads': sorted(counts)[len(counts) // 2] if counts else 0,
        'mean_reads': sum(counts) / len(counts) if counts else 0
    }

    stats['pct_ge_10'] = stats['junctions_ge_10'] / stats['total_junctions'] * 100 if counts else 0
    stats['pct_ge_20'] = stats['junctions_ge_20'] / stats['total_junctions'] * 100 if counts else 0

    return stats


def score_splice_sites(sequences_5ss, sequences_3ss):
    '''
    Score splice site strength using MaxEntScan.

    5' splice site (donor): 9bp (3 exon + 6 intron), typical score 8-10 bits
    3' splice site (acceptor): 23bp (20 intron + 3 exon), typical score 8-12 bits

    Requires: pip install maxentpy
    '''
    try:
        from maxentpy.maxent import score5, score3
    except ImportError:
        print('Install maxentpy: pip install maxentpy')
        return None, None

    scores_5ss = []
    for seq in sequences_5ss:
        if len(seq) == 9:
            try:
                scores_5ss.append(score5(seq.upper()))
            except Exception:
                scores_5ss.append(None)

    scores_3ss = []
    for seq in sequences_3ss:
        if len(seq) == 23:
            try:
                scores_3ss.append(score3(seq.upper()))
            except Exception:
                scores_3ss.append(None)

    return scores_5ss, scores_3ss


def generate_qc_report(bam_file, bed_file, output_prefix):
    '''
    Generate comprehensive splicing QC report.
    '''
    print(f'Analyzing: {bam_file}')
    print('=' * 50)

    # Junction saturation
    print('\n1. Junction Saturation Analysis')
    run_junction_saturation(bam_file, bed_file, output_prefix)

    # Junction annotation
    print('\n2. Junction Annotation')
    stats = run_junction_annotation(bam_file, bed_file, output_prefix)
    if stats is not None:
        print(stats)

    # Junction coverage
    print('\n3. Junction Read Coverage')
    junctions = count_junction_reads(bam_file)
    coverage_stats = analyze_junction_coverage(junctions)

    print(f"  Total junctions: {coverage_stats['total_junctions']}")
    print(f"  Junctions >= 10 reads: {coverage_stats['junctions_ge_10']} ({coverage_stats['pct_ge_10']:.1f}%)")
    print(f"  Junctions >= 20 reads: {coverage_stats['junctions_ge_20']} ({coverage_stats['pct_ge_20']:.1f}%)")
    print(f"  Median reads/junction: {coverage_stats['median_reads']}")

    # Quality assessment
    print('\n4. Quality Assessment')
    if coverage_stats['pct_ge_10'] > 50:
        print('  Junction coverage: GOOD (>50% junctions have >=10 reads)')
    elif coverage_stats['pct_ge_10'] > 30:
        print('  Junction coverage: ACCEPTABLE (30-50% junctions have >=10 reads)')
    else:
        print('  Junction coverage: POOR (<30% junctions have >=10 reads)')
        print('  Consider deeper sequencing for reliable splicing analysis')

    return coverage_stats


if __name__ == '__main__':
    # Example usage
    # generate_qc_report('sample.bam', 'annotation.bed', 'sample_qc')

    # Test splice site scoring
    # donor = 'CAGGTAAGT'  # Consensus 5'ss
    # acceptor = 'TTTTTTTTTTTTTTTTTTTTCAG'  # Example 3'ss
    # scores_5, scores_3 = score_splice_sites([donor], [acceptor])
    # print(f"5'ss score: {scores_5[0]:.2f}")
    # print(f"3'ss score: {scores_3[0]:.2f}")

    print('Provide BAM and BED files to run splicing QC analysis')
