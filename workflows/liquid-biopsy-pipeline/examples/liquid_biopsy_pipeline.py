#!/usr/bin/env python3
'''
Complete liquid biopsy analysis pipeline.
'''
# Reference: bwa 0.7.17+, vardict 1.8+, fgbio 2.1+, ichorcna 0.6.0+, numpy 1.26+, pandas 2.2+, pysam 0.22+, samtools 1.19+ | Verify API if version differs

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
    # bwa reads FASTQ; emit FASTQ carrying RX, align, then re-attach uBAM tags so RX survives for GroupReadsByUmi
    subprocess.run(f'samtools fastq -T RX {with_umis} | bwa mem -C -p -t 8 -Y {reference} - | '
                   f'fgbio ZipperBams --unmapped {with_umis} --ref {reference} --output {aligned}',
                   shell=True, check=True)

    sorted_bam = work_dir / f'{prefix}_sorted.bam'
    pysam.sort('-@', '8', '-o', str(sorted_bam), str(aligned))

    grouped = work_dir / f'{prefix}_grouped.bam'
    # --family-size-histogram is where genome-equivalents and the duplication plateau are read off.
    # A VAF without input GE is undefined (0.1% on 100 GE is noise; on 30,000 GE it is solid).
    subprocess.run([
        'fgbio', 'GroupReadsByUmi',
        '--input', str(sorted_bam),
        '--output', str(grouped),
        '--strategy', 'adjacency',
        '--edits', '1',
        '--family-size-histogram', str(work_dir / f'{prefix}_family_sizes.txt')
    ], check=True)

    consensus = work_dir / f'{prefix}_consensus.bam'
    subprocess.run([
        'fgbio', 'CallMolecularConsensusReads',
        '--input', str(grouped),
        '--output', str(consensus),
        '--min-reads', '1'
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


def load_panel_genes(bed_path):
    '''chrom -> [(start, end, gene)] from the panel BED. Column 4 is the gene, per VarDict's -g 4.'''
    regions = {}
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            f = line.rstrip('\n').split('\t')
            regions.setdefault(f[0], []).append((int(f[1]), int(f[2]), f[3] if len(f) > 3 else ''))
    return regions


def variants_to_df(vcf_path, panel_regions=None):
    '''Read a VCF into chrom/pos/ref/alt/gene rows.

    var2vcf_valid.pl declares no INFO/GENE tag and never writes the gene it parses, so reading
    rec.info['GENE'] silently yields '' for every record. Recover the gene from the same panel BED
    VarDict was handed via -g 4. BED is 0-based half-open, VCF pos is 1-based: start < pos <= end.
    '''
    rows = []
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for rec in vcf:
            gene = next((name for start, end, name in (panel_regions or {}).get(rec.chrom, [])
                         if start < rec.pos <= end), '')
            rows.append({'chrom': rec.chrom, 'pos': rec.pos, 'ref': rec.ref,
                         'alt': rec.alts[0] if rec.alts else '', 'gene': gene})
    return pd.DataFrame(rows)


def filter_chip(variants_df, wbc_variants=None):
    '''Subtract CHIP. Matched WBC is the definitive control; the gene list is a weak fallback.'''
    if wbc_variants is not None:
        # Any variant also seen in matched buffy-coat/WBC DNA is hematopoietic, not tumor -- regardless
        # of which gene it lands in. This is the commitment; the gene list below cannot substitute.
        key = ['chrom', 'pos', 'ref', 'alt']
        wbc_keys = set(map(tuple, wbc_variants[key].itertuples(index=False, name=None)))
        is_chip = variants_df[key].apply(lambda r: tuple(r) in wbc_keys, axis=1)
        return variants_df[~is_chip], variants_df[is_chip]

    # Fallback only. ~81.6% of cfDNA variants in controls and ~53.2% in cancer patients are CHIP
    # (Razavi 2019), and CHIP is not confined to these genes, so this UNDER-removes.
    print('WARNING: no matched WBC supplied -- gene-list CHIP filter only; calls are presumptively CHIP-contaminated')
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

        # CHIP subtraction is not optional: an unsubtracted plasma call is presumptively hematopoietic.
        panel = load_panel_genes(config['bed_file'])
        wbc = variants_to_df(config['wbc_vcf'], panel) if config.get('wbc_vcf') else None
        somatic, chip = filter_chip(variants_to_df(vcf, panel), wbc_variants=wbc)
        results['somatic'], results['chip'] = somatic, chip
        print(f'Somatic after CHIP subtraction: {len(somatic)}; CHIP removed: {len(chip)}')

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
    print()
    print("config keys: bam_file, reference, bed_file (col 4 = gene), output_vcf, data_type='panel'")
    print("             wbc_vcf  - matched buffy-coat/WBC calls; WITHOUT it, CHIP subtraction falls")
    print("                        back to a weak gene list and calls stay presumptively CHIP-contaminated")
