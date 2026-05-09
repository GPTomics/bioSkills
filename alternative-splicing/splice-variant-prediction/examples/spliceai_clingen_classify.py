#!/usr/bin/env python3
'''
SpliceAI VCF annotation with ClinGen SVI 2023 splicing-variant classification.

Workflow:
1. Run SpliceAI on a clinical VCF
2. Parse INFO field for delta scores (DS_AG, DS_AL, DS_DG, DS_DL)
3. Apply ClinGen SVI 2023 thresholds for PP3/BP4 evidence
4. Flag deep-intronic candidates for extended-window re-scoring
'''
# Reference: spliceai 1.3+, tensorflow 2.15+, pandas 2.2+ | Verify API if version differs

import re
import subprocess
import pandas as pd
from pathlib import Path


def run_spliceai(input_vcf, output_vcf, genome_fa, build='grch38', distance=50, mask=0):
    '''Run SpliceAI on a VCF.'''
    cmd = [
        'spliceai', '-I', str(input_vcf), '-O', str(output_vcf),
        '-R', str(genome_fa), '-A', build,
        '-D', str(distance), '-M', str(mask),
    ]
    subprocess.run(cmd, check=True)


def parse_spliceai_vcf(vcf_path):
    rows = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            chrom, pos, _, ref, alts, _, _, info = fields[:8]
            m = re.search(r'SpliceAI=([^;]+)', info)
            if not m:
                continue
            for ann in m.group(1).split(','):
                parts = ann.split('|')
                if len(parts) < 10:
                    continue
                allele, symbol = parts[0], parts[1]
                ds = [float(p) if p not in ('.', '') else 0.0 for p in parts[2:6]]
                rows.append({
                    'chrom': chrom, 'pos': int(pos), 'ref': ref, 'alt': allele,
                    'gene': symbol,
                    'DS_AG': ds[0], 'DS_AL': ds[1], 'DS_DG': ds[2], 'DS_DL': ds[3],
                    'delta_max': max(ds),
                })
    return pd.DataFrame(rows)


def apply_clingen_svi(df):
    '''Apply ClinGen SVI 2023 thresholds for PP3/BP4 (Walker 2023 AJHG).'''
    bins = [-0.01, 0.10, 0.20, 0.50, 0.80, 1.01]
    labels = ['BP4', 'inconclusive', 'PP3_supporting', 'PP3_moderate', 'PP3_strong']
    df = df.copy()
    df['acmg_evidence'] = pd.cut(df['delta_max'], bins=bins, labels=labels)
    return df


def flag_deep_intronic_candidates(df, threshold=0.05):
    '''Flag variants with weak default-window signal that might benefit from -D 2000 re-scoring.'''
    df = df.copy()
    df['extend_window_candidate'] = (df['delta_max'] >= threshold) & (df['delta_max'] < 0.20)
    return df


def main():
    input_vcf = Path('clinical_variants.vcf')
    spliceai_50 = Path('clinical_variants_spliceai_50.vcf')
    spliceai_2000 = Path('clinical_variants_spliceai_2000.vcf')
    genome = Path('GRCh38.primary_assembly.genome.fa')

    # Run default-window first for general clinical screening
    run_spliceai(input_vcf, spliceai_50, genome, distance=50)
    df_50 = parse_spliceai_vcf(spliceai_50)
    df_50 = apply_clingen_svi(df_50)
    df_50 = flag_deep_intronic_candidates(df_50)

    # Re-run with extended window only on candidates likely to be deep-intronic
    candidates = df_50[df_50['extend_window_candidate']]
    if not candidates.empty:
        candidates_vcf = Path('extend_window_candidates.vcf')
        candidates[['chrom', 'pos', 'ref', 'alt']].to_csv(candidates_vcf, sep='\t', index=False)
        run_spliceai(candidates_vcf, spliceai_2000, genome, distance=2000)

    df_50.to_csv('spliceai_clingen_classified.tsv', sep='\t', index=False)
    summary = df_50.groupby('acmg_evidence', observed=True).size()
    print(summary)


if __name__ == '__main__':
    main()
