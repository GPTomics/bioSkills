#!/usr/bin/env python3
'''
cfDNA preprocessing with UMI-aware deduplication using fgbio.
'''

import subprocess
import pysam
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def preprocess_cfdna_pipeline(input_bam, output_bam, reference, read_structure='3M2S+T', min_reads=2, threads=8):
    '''
    Full cfDNA preprocessing pipeline with fgbio.
    '''
    work_dir = Path(output_bam).parent
    prefix = Path(output_bam).stem

    with_umis = work_dir / f'{prefix}_umis.bam'
    subprocess.run([
        'fgbio', 'ExtractUmisFromBam',
        '--input', str(input_bam),
        '--output', str(with_umis),
        '--read-structure', read_structure, read_structure,
        '--single-tag', 'RX'
    ], check=True)

    aligned = work_dir / f'{prefix}_aligned.bam'
    cmd = f'bwa mem -t {threads} -Y {reference} {with_umis} | samtools view -bS - > {aligned}'
    subprocess.run(cmd, shell=True, check=True)

    sorted_bam = work_dir / f'{prefix}_sorted.bam'
    pysam.sort('-@', str(threads), '-o', str(sorted_bam), str(aligned))
    pysam.index(str(sorted_bam))

    grouped = work_dir / f'{prefix}_grouped.bam'
    subprocess.run([
        'fgbio', 'GroupReadsByUmi',
        '--input', str(sorted_bam),
        '--output', str(grouped),
        '--strategy', 'adjacency',
        '--edits', '1',
        '--min-map-q', '20'
    ], check=True)

    consensus = work_dir / f'{prefix}_consensus.bam'
    subprocess.run([
        'fgbio', 'CallMolecularConsensusReads',
        '--input', str(grouped),
        '--output', str(consensus),
        '--min-reads', str(min_reads),
        '--min-input-base-quality', '20'
    ], check=True)

    subprocess.run([
        'fgbio', 'FilterConsensusReads',
        '--input', str(consensus),
        '--output', str(output_bam),
        '--ref', str(reference),
        '--min-reads', str(min_reads),
        '--max-read-error-rate', '0.05'
    ], check=True)

    pysam.index(str(output_bam))
    return output_bam


def analyze_fragment_sizes(bam_path, max_size=500):
    '''Analyze cfDNA fragment size distribution.'''
    bam = pysam.AlignmentFile(bam_path, 'rb')
    sizes = []

    for read in bam.fetch():
        if read.is_proper_pair and not read.is_secondary and read.template_length > 0:
            if read.template_length <= max_size:
                sizes.append(read.template_length)

    bam.close()
    sizes = np.array(sizes)

    metrics = {
        'total_fragments': len(sizes),
        'median_size': np.median(sizes),
        'mean_size': np.mean(sizes),
        'modal_size': np.bincount(sizes).argmax() if len(sizes) > 0 else 0,
        'mononucleosome_fraction': np.sum((sizes >= 150) & (sizes <= 180)) / len(sizes) if len(sizes) > 0 else 0
    }

    return sizes, metrics


def plot_fragment_distribution(sizes, output_file):
    '''Plot cfDNA fragment size distribution.'''
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(sizes, bins=range(50, 400, 5), density=True, alpha=0.7, color='steelblue')
    ax.axvline(x=167, color='red', linestyle='--', label='Mononucleosome (~167bp)')
    ax.axvline(x=120, color='orange', linestyle=':', label='ctDNA enriched (~120bp)')

    ax.set_xlabel('Fragment Size (bp)')
    ax.set_ylabel('Density')
    ax.set_title('cfDNA Fragment Size Distribution')
    ax.legend()
    ax.set_xlim(50, 400)

    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    print('cfDNA Preprocessing Pipeline')
    print('=' * 40)
    print('1. preprocess_cfdna_pipeline() - Full UMI-aware pipeline')
    print('2. analyze_fragment_sizes() - Check fragment distribution')
    print('3. plot_fragment_distribution() - Visualize fragments')
