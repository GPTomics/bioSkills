# Reference: pandas 2.2+, scipy 1.12+, mageck 0.5.9+ | Verify API if version differs
#
# In vivo CRISPR screen per-animal meta-analysis.
# Combines per-animal MAGeCK output via Stouffer's Z method.

import pandas as pd
import numpy as np
from scipy.stats import norm
from pathlib import Path

# === INPUTS ===
# One mageck test gene_summary.txt per animal
animal_files = {
    'animal_1': 'animal_1.gene_summary.txt',
    'animal_2': 'animal_2.gene_summary.txt',
    'animal_3': 'animal_3.gene_summary.txt',
    'animal_4': 'animal_4.gene_summary.txt',
    'animal_5': 'animal_5.gene_summary.txt',
    'animal_6': 'animal_6.gene_summary.txt',
    'animal_7': 'animal_7.gene_summary.txt',
    'animal_8': 'animal_8.gene_summary.txt',
    'animal_9': 'animal_9.gene_summary.txt',
    'animal_10': 'animal_10.gene_summary.txt',
}

# === LOAD PER-ANIMAL RESULTS ===
animal_dfs = {}
for animal, file in animal_files.items():
    if Path(file).exists():
        animal_dfs[animal] = pd.read_csv(file, sep='\t')

if not animal_dfs:
    raise FileNotFoundError('No per-animal MAGeCK results found')

# === META-ANALYSIS VIA STOUFFER'S Z METHOD ===
# Z = sum(z_i) / sqrt(N); combines z-scores across animals
# (RRA p-values converted to z-scores)
all_results = pd.concat([df.assign(animal=name) for name, df in animal_dfs.items()])

def stouffer_meta(group):
    '''Stouffer Z method on negative-selection p-values.
    Returns: combined Z, p-value, animals significant.'''
    pvals = group['neg|p-value'].clip(lower=1e-10).values
    # Convert to z-scores; negative for dropout direction
    z_scores = -norm.ppf(pvals)  # one-sided
    combined_z = np.sum(z_scores) / np.sqrt(len(z_scores))
    combined_p = norm.sf(combined_z)
    return pd.Series({
        'meta_z': combined_z,
        'meta_p': combined_p,
        'mean_neg_score': group['neg|score'].mean(),
        'median_neg_lfc': group['neg|lfc'].median(),
        'n_animals': len(group),
        'animals_at_fdr_05': (group['neg|fdr'] < 0.05).sum(),
    })

meta = all_results.groupby('id').apply(stouffer_meta).reset_index()
meta = meta.sort_values('meta_z', ascending=False)

# === BH CORRECTION ON META P-VALUES ===
from statsmodels.stats.multitest import multipletests
meta['meta_fdr'] = multipletests(meta['meta_p'], method='fdr_bh')[1]

# === HIT CALLS ===
# Standard threshold: meta FDR <0.05 AND consistent across most animals (>50%)
n_animals = len(animal_dfs)
hits = meta[(meta['meta_fdr'] < 0.05) &
            (meta['animals_at_fdr_05'] >= n_animals * 0.5)]

print(f'Total animals analyzed: {n_animals}')
print(f'Meta-significant hits (FDR <0.05 + ≥50% animal consistency): {len(hits)}')
print(hits[['id', 'meta_z', 'meta_fdr', 'mean_neg_score',
             'animals_at_fdr_05', 'n_animals']].head(20).to_string(index=False))

# === EXPORT ===
hits.to_csv('in_vivo_meta_hits.tsv', sep='\t', index=False)
meta.to_csv('in_vivo_meta_all.tsv', sep='\t', index=False)

# === ALTERNATIVE: MAGeCK MLE WITH ANIMAL-AS-BATCH ===
print('\nAlternative analysis: MAGeCK MLE with animal-as-batch covariate')
print('  See SKILL.md for design matrix specification')
print('  This is preferred when running fewer animals (n=5-10) where Stouffer noisy')
