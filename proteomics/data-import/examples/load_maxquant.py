'''Load and clean a MaxQuant proteinGroups.txt: strip bookkeeping, resolve razor IDs, set 0->NaN, log2, diagnose missingness.

Self-contained: writes a tiny synthetic proteinGroups.txt to a tempdir and cleans up. No stray files.
'''
# Reference: pandas 2.2+, numpy 1.26+ | Verify API if version differs
import os
import tempfile
import pandas as pd
import numpy as np

def make_demo_protein_groups(path):
    # 10 real proteins spanning ~3 logs; below ~1e7 the protein falls under detection in S2,
    # so missingness is concentrated in the LOW-abundance proteins -- the left-censored MNAR signature.
    abund = [9.0e7, 6.0e7, 4.0e7, 2.5e7, 1.5e7, 9.0e6, 5.0e6, 2.0e6, 1.0e6, 6.0e5]
    rows = [{'Protein IDs': f'P{i:05d}', 'Gene names': f'GENE{i}', 'Reverse': '', 'Potential contaminant': '', 'Only identified by site': '',
             'LFQ intensity S1': a, 'LFQ intensity S2': (a * 0.95 if a >= 1.0e7 else 0), 'Intensity S1': a * 2, 'Intensity S2': a * 1.9}
            for i, a in enumerate(abund, 1)]
    rows += [
        {'Protein IDs': 'P12345;Q67890', 'Gene names': 'GENEA;GENEB', 'Reverse': '', 'Potential contaminant': '', 'Only identified by site': '', 'LFQ intensity S1': 1.2e7, 'LFQ intensity S2': 1.1e7, 'Intensity S1': 3.0e7, 'Intensity S2': 2.8e7},
        {'Protein IDs': 'REV__P99999', 'Gene names': '', 'Reverse': '+', 'Potential contaminant': '', 'Only identified by site': '', 'LFQ intensity S1': 2.0e6, 'LFQ intensity S2': 2.1e6, 'Intensity S1': 4.0e6, 'Intensity S2': 4.1e6},
        {'Protein IDs': 'CON__P02769', 'Gene names': 'ALB', 'Reverse': '', 'Potential contaminant': '+', 'Only identified by site': '', 'LFQ intensity S1': 8.0e7, 'LFQ intensity S2': 8.1e7, 'Intensity S1': 9.0e7, 'Intensity S2': 9.1e7},
        {'Protein IDs': 'P33333', 'Gene names': 'GENEC', 'Reverse': '', 'Potential contaminant': '', 'Only identified by site': '+', 'LFQ intensity S1': 1.0e6, 'LFQ intensity S2': 1.1e6, 'Intensity S1': 2.0e6, 'Intensity S2': 2.1e6},
    ]
    pd.DataFrame(rows).to_csv(path, sep='\t', index=False)

def load_clean_maxquant(path):
    pg = pd.read_csv(path, sep='\t', low_memory=False)  # mixed-type columns
    # All three flag columns are proteinGroups-only bookkeeping; '+' marks the row to drop
    mask = (pg.get('Reverse', '') != '+') & (pg.get('Potential contaminant', '') != '+') & (pg.get('Only identified by site', '') != '+')
    pg = pg[mask].copy()

    # Semicolon lists: first entry is the leading/razor identifier; Gene names may be blank
    pg['leading_protein'] = pg['Protein IDs'].str.split(';').str[0]
    pg['leading_gene'] = pg['Gene names'].where(pg['Gene names'].notna(), '').str.split(';').str[0]

    lfq_cols = [c for c in pg.columns if c.startswith('LFQ intensity ')]  # MaxLFQ-normalized, between-sample comparable
    matrix = pg[['leading_protein', 'leading_gene'] + lfq_cols].copy()
    matrix[lfq_cols] = matrix[lfq_cols].replace(0, np.nan)  # 0 means not-quantified; log2(0) = -inf
    matrix[lfq_cols] = np.log2(matrix[lfq_cols])
    return matrix, lfq_cols

def assess_missingness(matrix, sample_cols):
    total_pct = 100 * matrix[sample_cols].isna().sum().sum() / matrix[sample_cols].size
    mean_abund = matrix[sample_cols].mean(axis=1)  # negative corr with missingness => MNAR / left-censored
    mnar_corr = mean_abund.corr(matrix[sample_cols].isna().sum(axis=1))
    return total_pct, mnar_corr

with tempfile.TemporaryDirectory() as d:
    pg_path = os.path.join(d, 'proteinGroups.txt')
    make_demo_protein_groups(pg_path)
    matrix, lfq_cols = load_clean_maxquant(pg_path)
    print(f'Clean protein groups: {len(matrix)} (decoy/contaminant/site-only removed)')
    print(matrix.to_string(index=False))
    total_pct, mnar_corr = assess_missingness(matrix, lfq_cols)
    print(f'Missing: {total_pct:.1f}% | abundance-vs-missingness corr: {mnar_corr:.2f} (negative => MNAR)')
