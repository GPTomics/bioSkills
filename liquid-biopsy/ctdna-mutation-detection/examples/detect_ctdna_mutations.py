#!/usr/bin/env python3
'''
ctDNA mutation detection at low variant allele fractions.
'''

import subprocess
import pandas as pd
import pysam


CHIP_GENES = ['DNMT3A', 'TET2', 'ASXL1', 'PPM1D', 'JAK2', 'SF3B1', 'SRSF2', 'TP53', 'CBL', 'BCOR']


def call_variants_vardict(bam_file, reference, bed_file, output_vcf, min_vaf=0.005):
    '''Call variants with VarDict for high sensitivity.'''
    sample_id = bam_file.split('/')[-1].replace('.bam', '')

    cmd = f'''
    vardict-java \
        -G {reference} \
        -f {min_vaf} \
        -N {sample_id} \
        -b {bam_file} \
        -c 1 -S 2 -E 3 -g 4 \
        {bed_file} | \
    teststrandbias.R | \
    var2vcf_valid.pl \
        -N {sample_id} \
        -E \
        -f {min_vaf} \
        > {output_vcf}
    '''

    subprocess.run(cmd, shell=True, check=True)
    return output_vcf


def parse_vcf_variants(vcf_file):
    '''Parse VCF to DataFrame.'''
    variants = []

    with open(vcf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom, pos, _, ref, alt = parts[:5]
            info = dict(x.split('=') for x in parts[7].split(';') if '=' in x)

            variants.append({
                'chrom': chrom,
                'pos': int(pos),
                'ref': ref,
                'alt': alt,
                'af': float(info.get('AF', 0)),
                'dp': int(info.get('DP', 0)),
                'gene': info.get('GENE', '')
            })

    return pd.DataFrame(variants)


def filter_chip_variants(variants_df, chip_genes=None):
    '''Separate CHIP from somatic variants.'''
    if chip_genes is None:
        chip_genes = CHIP_GENES

    chip = variants_df[variants_df['gene'].isin(chip_genes)]
    somatic = variants_df[~variants_df['gene'].isin(chip_genes)]

    return somatic, chip


def track_mutations(bam_file, mutations, min_depth=100):
    '''Track specific mutations across samples.'''
    bam = pysam.AlignmentFile(bam_file, 'rb')
    results = []

    for chrom, pos, ref, alt in mutations:
        counts = {'ref': 0, 'alt': 0, 'other': 0}

        for pileupcolumn in bam.pileup(chrom, pos-1, pos):
            if pileupcolumn.pos != pos - 1:
                continue
            for read in pileupcolumn.pileups:
                if read.is_del or read.is_refskip:
                    continue
                base = read.alignment.query_sequence[read.query_position]
                if base == ref:
                    counts['ref'] += 1
                elif base == alt:
                    counts['alt'] += 1
                else:
                    counts['other'] += 1

        total = counts['ref'] + counts['alt'] + counts['other']
        vaf = counts['alt'] / total if total > 0 else 0

        results.append({
            'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
            'depth': total, 'alt_count': counts['alt'], 'vaf': vaf
        })

    bam.close()
    return pd.DataFrame(results)


if __name__ == '__main__':
    print('ctDNA Mutation Detection')
    print('=' * 40)
    print('1. call_variants_vardict() - High sensitivity calling')
    print('2. filter_chip_variants() - Remove CHIP mutations')
    print('3. track_mutations() - Track known variants')
    print()
    print('VAF detection limits:')
    print('  > 1%: Reliable with standard callers')
    print('  0.5-1%: Good with UMIs')
    print('  < 0.5%: Challenging, near noise floor')
