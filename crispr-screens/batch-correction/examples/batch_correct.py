import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

# === CONFIGURATION ===
output_dir = Path('batch_correction/')
output_dir.mkdir(exist_ok=True)

np.random.seed(42)

# === 1. SIMULATE MULTI-BATCH SCREEN DATA ===
print('Simulating multi-batch screen data...')

n_genes = 100
guides_per_gene = 4
n_guides = n_genes * guides_per_gene

gene_names = [f'Gene_{i:03d}' for i in range(n_genes)]
genes = np.repeat(gene_names, guides_per_gene)
guide_ids = [f'{g}_g{i+1}' for g in gene_names for i in range(guides_per_gene)]

batch1_baseline = 1000
batch2_baseline = 800

batch1_samples = ['B1_ctrl_1', 'B1_ctrl_2', 'B1_treat_1', 'B1_treat_2']
batch2_samples = ['B2_ctrl_1', 'B2_ctrl_2', 'B2_treat_1', 'B2_treat_2']
all_samples = batch1_samples + batch2_samples

counts = {}

essential_genes = gene_names[:10]
dropout_genes = gene_names[10:20]
hit_genes = gene_names[20:25]

for sample in batch1_samples:
    baseline = batch1_baseline + np.random.normal(0, 100)
    counts[sample] = np.random.poisson(baseline, n_guides).astype(float)

for sample in batch2_samples:
    baseline = batch2_baseline + np.random.normal(0, 100)
    counts[sample] = np.random.poisson(baseline, n_guides).astype(float)

for i, gene in enumerate(genes):
    if gene in essential_genes:
        for sample in all_samples:
            if 'treat' in sample:
                counts[sample][i] *= 0.3
    elif gene in hit_genes:
        for sample in all_samples:
            if 'treat' in sample:
                counts[sample][i] *= 0.4

counts_df = pd.DataFrame({'gene': genes, 'guide': guide_ids, **counts})
print(f'Simulated {n_guides} guides across {len(all_samples)} samples')

# === 2. ADD NON-TARGETING CONTROLS ===
n_controls = 50
for i in range(n_controls):
    control_row = {'gene': f'NonTargeting_{i+1:03d}', 'guide': f'NT_{i+1:03d}'}
    for sample in batch1_samples:
        control_row[sample] = np.random.poisson(batch1_baseline)
    for sample in batch2_samples:
        control_row[sample] = np.random.poisson(batch2_baseline)
    counts_df = pd.concat([counts_df, pd.DataFrame([control_row])], ignore_index=True)

print(f'Added {n_controls} non-targeting controls')

# === 3. MEDIAN NORMALIZATION ===
print('\n=== Median Normalization ===')

sample_cols = all_samples
sample_medians = counts_df[sample_cols].median()
global_median = sample_medians.median()
scale_factors = global_median / sample_medians

median_norm = counts_df.copy()
median_norm[sample_cols] = counts_df[sample_cols] * scale_factors

print('Scale factors:')
for s, f in scale_factors.items():
    print(f'  {s}: {f:.3f}')

# === 4. CONTROL-BASED NORMALIZATION ===
print('\n=== Control-Based Normalization ===')

nontargeting = counts_df['gene'].str.startswith('NonTargeting')
control_medians = counts_df.loc[nontargeting, sample_cols].median()
control_ref = control_medians.median()
control_factors = control_ref / control_medians

control_norm = counts_df.copy()
control_norm[sample_cols] = counts_df[sample_cols] * control_factors

print('Control-based scale factors:')
for s, f in control_factors.items():
    print(f'  {s}: {f:.3f}')

# === 5. CHECK BATCH EFFECTS ===
print('\n=== Batch Effect Analysis ===')

batch_labels = [1, 1, 1, 1, 2, 2, 2, 2]

from sklearn.decomposition import PCA

log_counts_raw = np.log2(counts_df[sample_cols].values.T + 1)
log_counts_norm = np.log2(control_norm[sample_cols].values.T + 1)

pca = PCA(n_components=2)
pcs_raw = pca.fit_transform(log_counts_raw)
var_explained_raw = pca.explained_variance_ratio_

pcs_norm = pca.fit_transform(log_counts_norm)
var_explained_norm = pca.explained_variance_ratio_

batch1_center_raw = pcs_raw[:4].mean(axis=0)
batch2_center_raw = pcs_raw[4:].mean(axis=0)
batch_dist_raw = np.linalg.norm(batch1_center_raw - batch2_center_raw)

batch1_center_norm = pcs_norm[:4].mean(axis=0)
batch2_center_norm = pcs_norm[4:].mean(axis=0)
batch_dist_norm = np.linalg.norm(batch1_center_norm - batch2_center_norm)

print(f'Batch separation (raw): {batch_dist_raw:.2f}')
print(f'Batch separation (normalized): {batch_dist_norm:.2f}')
print(f'Reduction: {(1 - batch_dist_norm/batch_dist_raw)*100:.1f}%')

# === 6. REPLICATE CORRELATION ===
print('\n=== Replicate Correlations ===')

replicate_pairs = [
    ('B1_ctrl_1', 'B1_ctrl_2'),
    ('B1_treat_1', 'B1_treat_2'),
    ('B2_ctrl_1', 'B2_ctrl_2'),
    ('B2_treat_1', 'B2_treat_2')
]

for r1, r2 in replicate_pairs:
    log_r1 = np.log2(control_norm[r1] + 1)
    log_r2 = np.log2(control_norm[r2] + 1)
    corr, pval = stats.pearsonr(log_r1, log_r2)
    print(f'  {r1} vs {r2}: r = {corr:.3f}')

# === 7. BATCH-AWARE LOG FOLD CHANGE ===
print('\n=== Batch-Aware Analysis ===')

def calculate_lfc(norm_df, treat_cols, ctrl_cols):
    treat_mean = norm_df[treat_cols].mean(axis=1)
    ctrl_mean = norm_df[ctrl_cols].mean(axis=1)
    return np.log2((treat_mean + 1) / (ctrl_mean + 1))

lfc_batch1 = calculate_lfc(control_norm, ['B1_treat_1', 'B1_treat_2'], ['B1_ctrl_1', 'B1_ctrl_2'])
lfc_batch2 = calculate_lfc(control_norm, ['B2_treat_1', 'B2_treat_2'], ['B2_ctrl_1', 'B2_ctrl_2'])

combined_lfc = (lfc_batch1 + lfc_batch2) / 2

control_norm['LFC_batch1'] = lfc_batch1
control_norm['LFC_batch2'] = lfc_batch2
control_norm['LFC_combined'] = combined_lfc

batch_corr, _ = stats.pearsonr(lfc_batch1, lfc_batch2)
print(f'LFC correlation between batches: {batch_corr:.3f}')

# === 8. SUMMARY STATISTICS ===
print('\n=== Hit Detection ===')

gene_lfc = control_norm.groupby('gene')['LFC_combined'].mean()

detected_hits = gene_lfc[gene_lfc < -1].index.tolist()
true_hits = essential_genes + list(hit_genes)

true_positives = len(set(detected_hits) & set(true_hits))
false_positives = len(set(detected_hits) - set(true_hits))

print(f'True hits: {len(true_hits)}')
print(f'Detected hits (LFC < -1): {len(detected_hits)}')
print(f'True positives: {true_positives}')
print(f'False positives: {false_positives}')

# === 9. EXPORT ===
print('\n=== Exporting Results ===')

counts_df.to_csv(output_dir / 'raw_counts.csv', index=False)
control_norm.to_csv(output_dir / 'normalized_counts.csv', index=False)

gene_summary = control_norm.groupby('gene').agg({
    'LFC_batch1': 'mean',
    'LFC_batch2': 'mean',
    'LFC_combined': 'mean'
}).reset_index()
gene_summary.to_csv(output_dir / 'gene_summary.csv', index=False)

qc_metrics = {
    'batch_separation_raw': batch_dist_raw,
    'batch_separation_norm': batch_dist_norm,
    'batch_lfc_correlation': batch_corr,
    'n_hits_detected': len(detected_hits),
    'true_positive_rate': true_positives / len(true_hits) if true_hits else 0
}
pd.DataFrame([qc_metrics]).to_csv(output_dir / 'qc_metrics.csv', index=False)

print(f'\nResults saved to {output_dir}/')
