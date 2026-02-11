'''Time-course analysis pipeline: temporal DE, clustering, optional periodicity, GAM fitting, enrichment.'''

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
from patsy import dmatrix
from sklearn.metrics import silhouette_score
from tslearn.clustering import TimeSeriesKMeans
from pygam import LinearGAM, s
import gseapy as gp

# --- Configuration ---
COUNTS_FILE = 'counts_normalized.csv'
METADATA_FILE = 'metadata.csv'
OUTPUT_PREFIX = 'timecourse_results'
# FDR <0.05: standard threshold; relax to 0.1 for exploratory clustering
FDR_THRESHOLD = 0.05
# 4-20 clusters typical; start with 8, evaluate with silhouette
N_CLUSTERS = 8
# softdtw gamma: lower = stricter alignment; 0.01 works for most time courses
DTW_GAMMA = 0.01
# GAM n_splines: 5 default; reduce to 3 for <6 time points
N_SPLINES = 5
# Set True if experiment covers 24h+ cycles with 2-4h resolution
CIRCADIAN_DESIGN = False

# --- Step 1: Load data ---
expr = pd.read_csv(COUNTS_FILE, index_col=0)
meta = pd.read_csv(METADATA_FILE)
print(f'Loaded: {expr.shape[0]} genes x {expr.shape[1]} samples, {meta["time"].nunique()} time points')

# --- Step 2: Temporal DE (F-test on spline basis) ---
# df=3: cubic spline; increase to 4-5 for >10 time points
spline_basis = dmatrix('bs(time, df=3)', data=meta, return_type='dataframe')
design_full = np.column_stack([np.ones(len(meta)), spline_basis.values])
design_reduced = np.ones((len(meta), 1))
df_diff = design_full.shape[1] - design_reduced.shape[1]
df_resid = len(meta) - design_full.shape[1]

pvals = []
for gene in expr.index:
    y = expr.loc[gene].values
    beta_full = np.linalg.lstsq(design_full, y, rcond=None)[0]
    beta_red = np.linalg.lstsq(design_reduced, y, rcond=None)[0]
    ss_full = np.sum((y - design_full @ beta_full) ** 2)
    ss_red = np.sum((y - design_reduced @ beta_red) ** 2)
    f_stat = ((ss_red - ss_full) / df_diff) / (ss_full / df_resid)
    pvals.append(1 - stats.f.cdf(f_stat, df_diff, df_resid))

_, fdr, _, _ = multipletests(pvals, method='fdr_bh')
temporal_mask = fdr < FDR_THRESHOLD
temporal_genes = expr.index[temporal_mask].tolist()
print(f'Significant temporal genes (FDR <{FDR_THRESHOLD}): {len(temporal_genes)}')

# QC gate: sufficient temporal genes
if len(temporal_genes) < 100:
    print('WARNING: Few temporal genes detected. Check replicates or relax FDR.')

expr_sig = expr.loc[temporal_genes]

# --- Step 3: Time-series clustering (tslearn soft-DTW) ---
# Standardize per-gene (row-wise): z-score each gene's profile across timepoints
expr_scaled = (expr_sig.values - expr_sig.values.mean(axis=1, keepdims=True)) / expr_sig.values.std(axis=1, keepdims=True)

# Reshape for tslearn: (n_genes, n_timepoints, 1)
X = expr_scaled.reshape(expr_scaled.shape[0], expr_scaled.shape[1], 1)

model = TimeSeriesKMeans(n_clusters=N_CLUSTERS, metric='softdtw',
                         metric_params={'gamma': DTW_GAMMA},
                         max_iter=50, random_state=42)
labels = model.fit_predict(X)

cluster_sizes = pd.Series(labels).value_counts().sort_index()
print('Cluster sizes:')
print(cluster_sizes.to_string())

# QC gate: no empty clusters
if (cluster_sizes == 0).any():
    print('WARNING: Empty clusters found. Reduce N_CLUSTERS.')

# Silhouette score for cluster quality validation
sil = silhouette_score(expr_scaled, labels, metric='euclidean')
print(f'Mean silhouette score: {sil:.3f}')

cluster_df = pd.DataFrame({'gene': temporal_genes, 'cluster': labels})
cluster_df.to_csv(f'{OUTPUT_PREFIX}_clusters.csv', index=False)

# --- Step 4a: Optional periodicity detection (CosinorPy) ---
if CIRCADIAN_DESIGN:
    from cosinorpy import cosinor

    # fit_group expects long-format DataFrame with columns 'x' (time), 'y' (expression), 'test' (gene name)
    # Map sample names to time values from metadata before reshaping
    # period=24: standard circadian; adjust for ultradian (4-12h) or infradian (>28h)
    records = []
    for gene in expr_sig.index:
        for sample_idx, sample in enumerate(expr_sig.columns):
            records.append({'x': meta['time'].iloc[sample_idx], 'y': expr_sig.loc[gene, sample], 'test': gene})
    expr_long = pd.DataFrame(records)
    cosinor_results = cosinor.fit_group(expr_long, period=24, n_components=1)
    rhythmic = cosinor_results[cosinor_results['p'] < 0.05]
    print(f'Rhythmic genes (p <0.05): {len(rhythmic)}')
    rhythmic.to_csv(f'{OUTPUT_PREFIX}_rhythmic.csv', index=False)

# --- Step 4b: GAM trajectory fitting ---
time_vals = meta['time'].values.reshape(-1, 1)

for cl_id in range(N_CLUSTERS):
    cl_mask = labels == cl_id
    if cl_mask.sum() == 0:
        continue

    mean_profile = expr_scaled[cl_mask].mean(axis=0)
    gam = LinearGAM(s(0, n_splines=N_SPLINES)).fit(time_vals, mean_profile)
    print(f'Cluster {cl_id}: GAM GCV = {gam.statistics_["GCV"]:.4f}')

# --- Step 5: Per-cluster pathway enrichment (gseapy/Enrichr) ---
clusters_with_terms = 0

for cl_id in range(N_CLUSTERS):
    cl_genes = cluster_df[cluster_df['cluster'] == cl_id]['gene'].tolist()
    if len(cl_genes) < 5:
        print(f'Cluster {cl_id}: too few genes ({len(cl_genes)}), skipping enrichment')
        continue

    enr = gp.enrichr(gene_list=cl_genes, gene_sets='GO_Biological_Process_2023',
                     organism='human', outdir=f'{OUTPUT_PREFIX}_enrichr_cluster_{cl_id}')
    sig_terms = enr.results[enr.results['Adjusted P-value'] < 0.05]
    n_terms = len(sig_terms)
    if n_terms > 0:
        clusters_with_terms += 1
    print(f'Cluster {cl_id}: {n_terms} significant GO terms')

# QC gate: enrichment coverage
print(f'Clusters with significant GO terms: {clusters_with_terms} / {N_CLUSTERS}')
if clusters_with_terms < 3:
    print('WARNING: Few clusters enriched. Check gene naming convention or relax thresholds.')

print(f'\nPipeline complete: {len(temporal_genes)} temporal genes, {N_CLUSTERS} clusters, {clusters_with_terms} enriched')
