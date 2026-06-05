#!/usr/bin/env python3
'''
ctDNA mutation detection: de novo VarDict calling, known-variant VAF tracking, and CHIP-gene flagging.

Detection (tracking a known variant set) and de novo calling are different statistical problems.
VarDict/VEP are external CLIs verified separately against their docs; this script wraps the call
and implements the pysam-based tracker and CHIP filter in Python.
'''
# Reference: pysam 0.22+, pandas 2.2+ | Verify API if version differs

import subprocess
import pandas as pd
import pysam

# Canonical CHIP genes (Razavi 2019 Nat Med 25:1928; Jaiswal 2014 NEJM 371:2488).
# A flag for extra scrutiny only -- NOT a substitute for matched-WBC subtraction,
# since TP53/ATM are also bona fide tumor suppressors.
CHIP_GENES = ['DNMT3A', 'TET2', 'ASXL1', 'PPM1D', 'TP53', 'JAK2', 'SF3B1', 'SRSF2', 'GNB1', 'GNAS', 'CBL', 'ATM', 'CHEK2']


def call_variants_vardict(bam_file, reference, bed_file, output_vcf, min_vaf=0.005):
    '''Run de novo low-VAF calling on a UMI-consensus BAM.

    min_vaf 0.005 = 0.5%, the practical UMI-consensus de novo floor; below this the alt
    signal approaches the per-base error floor. var2vcf_valid.pl -f must match VarDict -f
    or calls are silently re-filtered. -c/-S/-E/-g are 1-based BED column indices, not
    coordinates. var2vcf_valid.pl -E suppresses the END tag (opposite of VarDict's -E).
    '''
    sample_id = bam_file.split('/')[-1].replace('.bam', '')
    cmd = (
        f'vardict-java -G {reference} -f {min_vaf} -N {sample_id} -b {bam_file} '
        f'-c 1 -S 2 -E 3 -g 4 {bed_file} | '
        f'teststrandbias.R | '
        f'var2vcf_valid.pl -N {sample_id} -E -f {min_vaf} > {output_vcf}'
    )
    subprocess.run(cmd, shell=True, check=True)
    return output_vcf


def parse_vcf_variants(vcf_file):
    '''Parse a VarDict VCF into a DataFrame of chrom/pos/ref/alt/af/dp/gene.'''
    rows = []
    with open(vcf_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            chrom, pos, _, ref, alt = parts[:5]
            info = dict(kv.split('=') for kv in parts[7].split(';') if '=' in kv)
            rows.append({'chrom': chrom, 'pos': int(pos), 'ref': ref, 'alt': alt,
                         'af': float(info.get('AF', 0)), 'dp': int(info.get('DP', 0)),
                         'gene': info.get('GENE', '')})
    return pd.DataFrame(rows)


def flag_chip_variants(variants_df, chip_genes=None):
    '''Split calls into candidate-somatic and CHIP-gene subsets.

    Heuristic gene flag only. Confident somatic-tumor classification requires matched-WBC
    subtraction (a cfDNA variant present in WBC is CHIP or germline -> remove).
    '''
    genes = chip_genes if chip_genes is not None else CHIP_GENES
    is_chip = variants_df['gene'].isin(genes)
    return variants_df[~is_chip], variants_df[is_chip]


def track_known_variants(bam_file, variants):
    '''Quantify VAF of a pre-specified variant set at fixed loci (the detection problem).

    variants: list of (chrom, pos, ref, alt) with pos 1-based. truncate=True keeps only the
    requested column; aggregate alt_count across loci as the panel-level detection signal.
    SNV loci only - the single-base compare scores every indel read as other and would
    report an indel reporter as cleared; indels need read.indel/CIGAR-aware counting.
    '''
    bam = pysam.AlignmentFile(bam_file, 'rb')
    rows = []
    for chrom, pos, ref, alt in variants:
        counts = {'ref': 0, 'alt': 0, 'other': 0}
        for col in bam.pileup(chrom, pos - 1, pos, truncate=True):
            for read in col.pileups:
                if read.is_del or read.is_refskip:
                    continue
                base = read.alignment.query_sequence[read.query_position]
                counts['alt' if base == alt else 'ref' if base == ref else 'other'] += 1
        depth = counts['ref'] + counts['alt'] + counts['other']
        rows.append({'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
                     'depth': depth, 'alt_count': counts['alt'],
                     'vaf': counts['alt'] / depth if depth else 0.0})
    bam.close()
    return pd.DataFrame(rows)


if __name__ == '__main__':
    print('ctDNA mutation detection')
    print('=' * 40)
    print('call_variants_vardict() - de novo low-VAF calling on a consensus BAM')
    print('parse_vcf_variants()    - VarDict VCF -> DataFrame')
    print('flag_chip_variants()    - heuristic CHIP-gene split (matched WBC still required)')
    print('track_known_variants()  - pileup VAF at a pre-specified locus set (MRD)')
    print()
    print('VAF regime -> approach:')
    print('  > 1%      : any caller on a deduplicated BAM')
    print('  0.1-1%    : single-strand UMI consensus + VarDict/Mutect2/umi-varcal')
    print('  < 0.1%    : duplex consensus + tumor-informed integration')
