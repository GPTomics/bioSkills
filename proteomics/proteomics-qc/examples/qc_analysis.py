'''Three-level proteomics QC: inspect raw signal and remove contaminants BEFORE normalizing.

Self-contained: synthesizes a small protein matrix with one deliberately low-loaded
sample and one contaminant row, then runs experiment-level QC the way the skill mandates.
Writes nothing to disk.'''
# Reference: pandas 2.2+, numpy 1.26+, scikit-learn 1.4+ | Verify API if version differs
import numpy as np
import pandas as pd
from itertools import combinations
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.stats import f_oneway

rng = np.random.default_rng(0)

n_proteins = 400
samples = [f'{cond}_{rep}' for cond in ['ctrl', 'treat'] for rep in range(1, 4)]
sample_groups = pd.Series([s.split('_')[0] for s in samples], index=samples)

# log-normal abundances; treat shifts a subset of proteins
base = rng.normal(20, 2, n_proteins)
linear = pd.DataFrame({s: 2.0 ** (base + rng.normal(0, 0.3, n_proteins)) for s in samples},
                      index=[f'P{i:04d}' for i in range(n_proteins)])
treat_up = rng.choice(n_proteins, 40, replace=False)
for s in [c for c in samples if c.startswith('treat')]:
    linear.iloc[treat_up, linear.columns.get_loc(s)] *= 3.0

# inject realistic faults: ctrl_3 loaded ~3x low; one dominant keratin contaminant row
linear['ctrl_3'] *= 0.33
linear.loc['CON_KRT1'] = linear.iloc[0].values * 40.0

# MNAR-style missingness: drop the lowest-abundance values
threshold = linear.stack().quantile(0.10)
linear = linear.mask(linear < threshold)

protein_groups = linear.copy()
protein_groups['Potential contaminant'] = ['+' if idx.startswith('CON_') else '' for idx in protein_groups.index]
protein_groups['Reverse'] = ''
intensity_cols = samples

print('=== Level 3 QC: inspect-before-normalize ===\n')

# Step 1: RAW per-sample signal BEFORE any transform -- this is where loading failures are visible
raw = protein_groups[intensity_cols]
raw_qc = pd.DataFrame({'n_quantified': raw.notna().sum(),
                       'total_signal': raw.sum(),
                       'median_intensity': raw.median()})
group_median_total = raw_qc['total_signal'].median()
raw_qc['fold_vs_median'] = raw_qc['total_signal'] / group_median_total
print('Raw per-sample signal (loading failures show HERE, before normalization):')
print(raw_qc.round(2), '\n')
loading_failures = raw_qc.index[raw_qc['fold_vs_median'] < 0.5].tolist()  # <0.5x group median = loading/injection failure
print(f'Flagged loading failures (>2x low): {loading_failures}\n')

# Step 2: contaminant fraction, then strip contaminant/decoy rows BEFORE log + normalize
contaminant_flags = ['Potential contaminant', 'Reverse', 'Only identified by site']
keep = pd.Series(True, index=protein_groups.index)
for col in contaminant_flags:
    match = next((c for c in protein_groups.columns if c.lower() == col.lower()), None)
    if match is not None:
        keep &= protein_groups[match].fillna('') != '+'
contaminant_frac = 100 * raw[~keep.values].sum().sum() / raw.sum().sum()
print(f'Contaminant fraction of summed intensity: {contaminant_frac:.1f}% (PTXQC flags >1%)')
clean = protein_groups[keep][intensity_cols]
print(f'Rows after contaminant removal: {len(clean)} (was {len(protein_groups)})\n')

# Step 3: NOW log2-transform and median-normalize on survivors (mechanics route to quantification)
log2 = np.log2(clean)
survivors = [s for s in intensity_cols if s not in loading_failures]
normalized = log2[survivors] - log2[survivors].median()

# Replicate correlation on log2
corr = normalized.corr(method='pearson')
print('Within-group replicate correlation (log2):')
groups_kept = sample_groups[survivors]
for group in groups_kept.unique():
    members = groups_kept[groups_kept == group].index
    for s1, s2 in combinations(members, 2):
        print(f'  {group}: r({s1}, {s2}) = {corr.loc[s1, s2]:.3f}')

# CV on the LINEAR scale (base formula is only valid on linear intensity)
print('\nMedian CV per condition (linear scale):')
for group in groups_kept.unique():
    block = clean[groups_kept[groups_kept == group].index]
    cv = block.std(axis=1) / block.mean(axis=1)
    print(f'  {group}: {100 * cv.median():.1f}%')

# Missingness mechanism: present-fraction vs abundance (rising = MNAR)
present_frac = normalized.notna().mean(axis=1)
mean_abundance = normalized.mean(axis=1)
low, high = mean_abundance.quantile(0.25), mean_abundance.quantile(0.75)
print(f'\nMissingness: present-fraction low-abundance={present_frac[mean_abundance <= low].mean():.2f} '
      f'high-abundance={present_frac[mean_abundance >= high].mean():.2f} (rising-with-abundance = MNAR -> impute LOW)')

# PCA / batch (impute temporarily for projection only; drop rows missing across all survivors)
complete = normalized.dropna(how='all')
imputed = complete.apply(lambda r: r.fillna(r.median()), axis=1)
scaled = StandardScaler().fit_transform(imputed.T)
pca = PCA(n_components=3).fit(scaled)
coords = pd.DataFrame(pca.transform(scaled), columns=['PC1', 'PC2', 'PC3'], index=survivors)
coords['condition'] = groups_kept.values
print('\nPCA variance explained:', [f'{100 * v:.1f}%' for v in pca.explained_variance_ratio_])
cond_groups = [coords[coords['condition'] == c]['PC1'] for c in coords['condition'].unique()]
_, p = f_oneway(*cond_groups)
print(f'PC1 ~ condition: p={p:.4f} (low p = biology drives PC1, good)')

print('\n=== Summary ===')
print(f'Excluded loading failures: {loading_failures}')
print(f'Removed {len(protein_groups) - len(clean)} contaminant/decoy rows before normalization')
print('Normalized and ran correlation/CV/PCA on survivors only')
