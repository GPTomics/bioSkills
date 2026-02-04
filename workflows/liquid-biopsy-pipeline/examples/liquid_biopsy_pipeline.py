#!/usr/bin/env python3
'''
Complete liquid biopsy analysis pipeline.
'''

import subprocess
import pysam
import numpy as np
import pandas as pd
from pathlib import Path


def check_preanalytical_quality(sample_metadata):
    '''Check pre-analytical factors.'''
    issues = []
    if sample_metadata.get('tube_type') == 'EDTA':
        if sample_metadata.get('processing_delay_hours', 0) > 6:
            issues.append('EDTA processed > 6 hours')
    if sample_metadata.get('hemolysis_score', 0) > 1:
        issues.append('Hemolysis detected')
    return issues


def preprocess_cfdna(input_bam, output_bam, reference, work_dir):
    '''UMI-aware preprocessing with fgbio.'''
    work_dir = Path(work_dir)
    prefix = Path(output_bam).stem

    with_umis = work_dir / f'{prefix}_umis.bam'
    subprocess.run([
        'fgbio', 'ExtractUmisFromBam',
        '--input', str(input_bam),
        '--output', str(with_umis),
        '--read-structure', '3M2S+T', '3M2S+T',
        '--single-tag', 'RX'
    ], check=True)

    aligned = work_dir / f'{prefix}_aligned.bam'
    subprocess.run(f'bwa mem -t 8 -Y {reference} {with_umis} | samtools view -bS - > {aligned}',
                   shell=True, check=True)

    sorted_bam = work_dir / f'{prefix}_sorted.bam'
    pysam.sort('-@', '8', '-o', str(sorted_bam), str(aligned))

    grouped = work_dir / f'{prefix}_grouped.bam'
    subprocess.run([
        'fgbio', 'GroupReadsByUmi',
        '--input', str(sorted_bam),
        '--output', str(grouped),
        '--strategy', 'adjacency'
    ], check=True)

    consensus = work_dir / f'{prefix}_consensus.bam'
    subprocess.run([
        'fgbio', 'CallMolecularConsensusReads',
        '--input', str(grouped),
        '--output', str(consensus),
        '--min-reads', '2'
    ], check=True)

    subprocess.run([
        'fgbio', 'FilterConsensusReads',
        '--input', str(consensus),
        '--output', str(output_bam),
        '--ref', str(reference),
        '--min-reads', '2'
    ], check=True)

    return output_bam


def verify_cfdna_quality(bam_path):
    '''Verify cfDNA fragment profile.'''
    bam = pysam.AlignmentFile(bam_path, 'rb')
    sizes = []
    for read in bam.fetch():
        if read.is_proper_pair and not read.is_secondary and 0 < read.template_length <= 400:
            sizes.append(read.template_length)
    bam.close()

    sizes = np.array(sizes)
    modal = np.bincount(sizes).argmax() if len(sizes) > 0 else 0
    mono_frac = np.sum((sizes >= 150) & (sizes <= 180)) / len(sizes) if len(sizes) > 0 else 0
    qc_pass = 150 <= modal <= 180 and mono_frac > 0.3

    return {'modal_size': modal, 'mono_fraction': mono_frac, 'qc_pass': qc_pass}


def call_variants_vardict(bam_file, reference, bed_file, output_vcf, min_vaf=0.005):
    '''Call variants with VarDict.'''
    sample_id = Path(bam_file).stem
    cmd = f'''
    vardict-java -G {reference} -f {min_vaf} -N {sample_id} -b {bam_file} \
        -c 1 -S 2 -E 3 -g 4 {bed_file} | \
    teststrandbias.R | var2vcf_valid.pl -N {sample_id} -E -f {min_vaf} > {output_vcf}
    '''
    subprocess.run(cmd, shell=True, check=True)
    return output_vcf


def filter_chip(variants_df):
    '''Filter CHIP variants.'''
    chip_genes = ['DNMT3A', 'TET2', 'ASXL1', 'PPM1D', 'JAK2', 'SF3B1', 'SRSF2', 'TP53']
    chip = variants_df[variants_df['gene'].isin(chip_genes)]
    somatic = variants_df[~variants_df['gene'].isin(chip_genes)]
    return somatic, chip


def run_pipeline(config):
    '''Run complete liquid biopsy pipeline.'''
    results = {}

    # Check pre-analytical
    if 'metadata' in config:
        issues = check_preanalytical_quality(config['metadata'])
        if issues:
            print(f'Pre-analytical issues: {issues}')
        results['preanalytical_issues'] = issues

    # Preprocess
    if config.get('has_umis'):
        bam = preprocess_cfdna(config['bam_file'], config['output_bam'],
                               config['reference'], config['work_dir'])
    else:
        bam = config['bam_file']

    # Fragment QC
    frag_qc = verify_cfdna_quality(bam)
    results['fragment_qc'] = frag_qc
    if not frag_qc['qc_pass']:
        print(f"WARNING: Atypical fragment profile (modal: {frag_qc['modal_size']}bp)")

    # Analysis based on data type
    if config['data_type'] == 'panel':
        vcf = call_variants_vardict(bam, config['reference'], config['bed_file'],
                                    config['output_vcf'])
        results['vcf'] = vcf

    return results


if __name__ == '__main__':
    print('Liquid Biopsy Pipeline')
    print('=' * 40)
    print('1. check_preanalytical_quality() - Pre-analytical QC')
    print('2. preprocess_cfdna() - UMI preprocessing')
    print('3. verify_cfdna_quality() - Fragment QC')
    print('4. call_variants_vardict() - Mutation detection')
    print('5. filter_chip() - Remove CHIP variants')
    print('6. run_pipeline() - Complete pipeline')
