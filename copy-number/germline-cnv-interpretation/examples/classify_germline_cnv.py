#!/usr/bin/env python3
'''Triage ClassifyCNV output against the 2019 ACMG/ClinGen points framework.

ClassifyCNV scores the automatable evidence (gene content, ClinGen dosage overlap,
population frequency). It cannot score de novo status, segregation, or literature --
so unsupervised output systematically lands in VUS. This script flags the VUS CNVs
sitting near a tier boundary, where adding case-specific Section 4-5 evidence would
change the classification.

Run ClassifyCNV first:
  python ClassifyCNV.py --infile cnvs.bed --GenomeBuild hg38 --precise \
      --outdir classifycnv_out
'''
# Reference: ClassifyCNV 1.1+, pandas 2.2+ | Verify API if version differs

import sys
import pandas as pd

# 2019 ACMG/ClinGen tier boundaries (Riggs 2020). Total score maps to classification.
TIERS = [
    (0.99, float('inf'), 'Pathogenic'),
    (0.90, 0.98, 'Likely pathogenic'),
    (-0.89, 0.89, 'VUS'),
    (-0.98, -0.90, 'Likely benign'),
    (float('-inf'), -0.99, 'Benign'),
]


def classify_from_score(score):
    for lo, hi, label in TIERS:
        if lo <= score <= hi:
            return label
    return 'VUS'


def triage(scoresheet, out='cnv_triage.tsv'):
    '''Flag VUS CNVs near a tier boundary that need manual Section 4-5 evidence.'''
    df = pd.read_csv(scoresheet, sep='\t')
    score_col = 'Total score' if 'Total score' in df.columns else 'Score'

    # A VUS within ~0.3 points of either likely-pathogenic or likely-benign is where
    # de novo / segregation evidence (not scored by the tool) is decisive.
    df['toward_pathogenic'] = df[score_col].between(0.60, 0.89)
    df['toward_benign'] = df[score_col].between(-0.89, -0.60)
    df['needs_manual_evidence'] = (
        (df.get('Classification', df[score_col].map(classify_from_score)) == 'VUS')
        & (df['toward_pathogenic'] | df['toward_benign']))

    n_flag = int(df['needs_manual_evidence'].sum())
    df.to_csv(out, sep='\t', index=False)
    print(f'{len(df)} CNVs scored; {n_flag} VUS near a tier boundary need manual '
          f'de novo / segregation / literature evidence before a final classification.')
    print(f'Output: {out}')
    return df


if __name__ == '__main__':
    sheet = sys.argv[1] if len(sys.argv) > 1 else 'classifycnv_out/Scoresheet.txt'
    triage(sheet)
