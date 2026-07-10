'''Time-course pipeline: temporal DE -> soft clustering -> GAM trajectory -> per-cluster enrichment.

Self-contained demo: generates a small synthetic time course, runs the pipeline end-to-end,
and writes all outputs to a temporary directory that is removed on exit (leaves no stray files).
The circadian rhythm branch is OPTIONAL and gated OFF by default (see CIRCADIAN_DESIGN).
'''
# Reference: numpy 1.26+, pandas 2.2+, scipy 1.12+, statsmodels 0.14+, patsy 1.0+, scikit-learn 1.4+, tslearn 0.6+, pygam 0.9+, gseapy 1.2+ | Verify API if version differs

import shutil
import tempfile
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from patsy import dmatrix
from sklearn.metrics import silhouette_score
from tslearn.clustering import TimeSeriesKMeans
from pygam import LinearGAM, s
import gseapy as gp

# --- Configuration ---
FDR_THRESHOLD = 0.05      # standard temporal-DE threshold; 0.1 for exploratory clustering only
N_CLUSTERS = 4            # a CHOICE not a result; matches the 4 synthetic archetypes here, sweep on real data
DTW_GAMMA = 0.01          # soft-DTW smoothing; lower = closer to hard DTW alignment
N_SPLINES = 5             # GAM basis-dimension ceiling; keep < number of unique timepoints
CIRCADIAN_DESIGN = False  # set True ONLY under a real circadian design (see the gate below)
RANDOM_STATE = 42

# --- Step 0: synthetic data (8 timepoints x 3 replicates; 4 temporal archetypes + flat genes) ---
rng = np.random.RandomState(RANDOM_STATE)
timepoints = np.array([0, 3, 6, 9, 12, 18, 24, 36])
n_reps = 3
times = np.repeat(timepoints, n_reps)
tnorm = (timepoints - timepoints.mean()) / timepoints.std()
archetypes = {'early_transient': np.exp(-((tnorm + 1) ** 2)),
              'late_sustained': 1 / (1 + np.exp(-3 * tnorm)),
              'monotone_down': -tnorm,
              'biphasic': np.sin(1.5 * tnorm)}
rows, gene_ids, gene_program = {}, [], {}
for prog, shape in archetypes.items():
    for k in range(40):
        gid = f'{prog}_{k}'
        profile = np.repeat(2.0 * shape, n_reps) + 0.3 * rng.randn(len(times))
        rows[gid] = profile
        gene_ids.append(gid)
        gene_program[gid] = prog
for k in range(80):  # flat / non-temporal genes: pure noise around a random level
    gid = f'flat_{k}'
    rows[gid] = rng.uniform(3, 7) + 0.3 * rng.randn(len(times))
    gene_ids.append(gid)
expr = pd.DataFrame(rows).T
expr.columns = [f's{i}' for i in range(len(times))]
meta = pd.DataFrame({'sample': expr.columns, 'time': times})
print(f'Synthetic input: {expr.shape[0]} genes x {expr.shape[1]} samples, {meta["time"].nunique()} timepoints')

workdir = tempfile.mkdtemp(prefix='timecourse_')

# --- Step 1: temporal DE (F-test of a spline model vs intercept) ---
spline_basis = dmatrix('bs(time, df=3)', data=meta, return_type='dataframe')  # df=3 cubic; 4-5 for >10 timepoints
# patsy already includes an Intercept column; do NOT prepend another np.ones (duplicate intercept
# inflates model df and miscalibrates the F-test -> wrong temporal-DE FDR). Use the basis directly.
design_full = spline_basis.values                          # [Intercept, bs1, bs2, bs3] = 4 cols
design_reduced = np.ones((len(meta), 1))
df_diff = design_full.shape[1] - design_reduced.shape[1]   # 4 - 1 = 3
df_resid = len(meta) - design_full.shape[1]                # n - 4

pvals = []
for gene in expr.index:
    y = expr.loc[gene].values
    ss_full = np.sum((y - design_full @ np.linalg.lstsq(design_full, y, rcond=None)[0]) ** 2)
    ss_red = np.sum((y - design_reduced @ np.linalg.lstsq(design_reduced, y, rcond=None)[0]) ** 2)
    f_stat = ((ss_red - ss_full) / df_diff) / (ss_full / df_resid)
    pvals.append(1 - stats.f.cdf(f_stat, df_diff, df_resid))

# multipletests default is Holm-Sidak; force BH explicitly
_, fdr, _, _ = multipletests(pvals, method='fdr_bh')
temporal_genes = expr.index[fdr < FDR_THRESHOLD].tolist()
print(f'Significant temporal genes (FDR <{FDR_THRESHOLD}): {len(temporal_genes)}')
if len(temporal_genes) < 100:
    print('WARNING: Few temporal genes. On real data, check replicates or relax FDR.')

# --- Step 2: filter to the temporal genes (clustering input; never the full matrix) ---
expr_sig = expr.loc[temporal_genes]

# --- Step 3: soft clustering on per-gene z-scored profiles (collapse replicates to timepoint means) ---
sig_means = expr_sig.T.groupby(meta['time'].values).mean().T   # gene x timepoint means
scaled = (sig_means.values - sig_means.values.mean(axis=1, keepdims=True)) / sig_means.values.std(axis=1, keepdims=True)
X = scaled.reshape(scaled.shape[0], scaled.shape[1], 1)

# soft-DTW tolerates phase-shifted profiles; use 'euclidean' when absolute phase is meaningful
model = TimeSeriesKMeans(n_clusters=N_CLUSTERS, metric='softdtw', metric_params={'gamma': DTW_GAMMA},
                         max_iter=50, random_state=RANDOM_STATE)
labels = model.fit_predict(X)
cluster_df = pd.DataFrame({'gene': temporal_genes, 'cluster': labels})
cluster_df.to_csv(f'{workdir}/clusters.csv', index=False)

sizes = pd.Series(labels).value_counts().sort_index()
print('Cluster sizes:\n' + sizes.to_string())
if (sizes == 0).any():
    print('WARNING: Empty clusters. Reduce N_CLUSTERS.')
# silhouette on the same (Euclidean, z-scored) geometry used to summarize shape
print(f'Mean silhouette: {silhouette_score(scaled, labels, metric="euclidean"):.3f}')

# --- Step 4a: OPTIONAL rhythm detection - GATED (skipped unless the design licenses it) ---
n_cycles = (timepoints.max() - timepoints.min()) / 24.0
samples_per_cycle = len(timepoints) / max(n_cycles, 1e-9)
# Two independent preconditions: the COMPUTABLE gate (cycles/samples) and the un-computable affirmation
# (CIRCADIAN_DESIGN) that the run is a circadian design with randomized collection order. Report which failed.
gate_design = n_cycles >= 2 and samples_per_cycle >= 6
if CIRCADIAN_DESIGN and gate_design:
    from CosinorPy import cosinor   # import name is capitalized CosinorPy, not cosinorpy
    records = [{'x': meta['time'].iloc[j], 'y': expr_sig.iloc[i, j], 'test': g}
               for i, g in enumerate(expr_sig.index) for j in range(expr_sig.shape[1])]
    res = cosinor.fit_group(pd.DataFrame(records), period=24, n_components=1)
    # fit_group returns a BH-adjusted 'q' column across genes; use it (not raw 'p') plus a RELATIVE-amplitude
    # filter (rAMP = amplitude/mesor; fit_group has no rAMP column) -- significance alone over-detects rhythms.
    res['rAMP'] = res['amplitude'] / res['mesor']
    rhythmic = res[(res['q'] < 0.05) & (res['rAMP'] > 0.1)]
    print(f'Rhythmic genes (BH q <0.05 and rAMP >0.1): {len(rhythmic)}')
elif not gate_design:
    print(f'Rhythm detection SKIPPED: sampling design is inadequate - {n_cycles:.1f} cycles, '
          f'{samples_per_cycle:.1f} samples/cycle (need >=2 cycles and >=6-8/cycle).')
else:
    print('Rhythm detection SKIPPED: design meets the cycle/sampling floor but CIRCADIAN_DESIGN is not set. '
          'Set it True ONLY for a circadian design with randomized collection order (the un-computable precondition).')

# --- Step 4b: GAM trajectory per cluster (on standardized cluster-mean profiles -> Gaussian is fine) ---
tvals = np.unique(times).reshape(-1, 1)
for cl_id in range(N_CLUSTERS):
    if (labels == cl_id).sum() == 0:
        continue
    mean_profile = scaled[labels == cl_id].mean(axis=0)
    # n_splines is the basis-dimension CEILING (like mgcv k); the penalty picks realized wiggliness (edof)
    gam = LinearGAM(s(0, n_splines=N_SPLINES)).fit(tvals, mean_profile)
    print(f'Cluster {cl_id}: GAM GCV = {gam.statistics_["GCV"]:.4f}, edof = {gam.statistics_["edof"]:.2f}')

# --- Step 5: per-cluster enrichment (background = temporal genes, NOT the genome) ---
# gp.enrich() runs a local hypergeometric test with an explicit background (offline, no files).
# On real data with gene symbols, gp.enrichr(..., background=temporal_genes, outdir=None) queries GO libraries.
gene_sets = {prog: [g for g in temporal_genes if gene_program.get(g) == prog] for prog in archetypes}
background = list(expr_sig.index)
clusters_with_terms = 0
for cl_id in range(N_CLUSTERS):
    cl_genes = cluster_df[cluster_df['cluster'] == cl_id]['gene'].tolist()
    if len(cl_genes) < 5:
        print(f'Cluster {cl_id}: too few genes ({len(cl_genes)}), skipping enrichment')
        continue
    enr = gp.enrich(gene_list=cl_genes, gene_sets=gene_sets, background=background, outdir=None, no_plot=True)
    sig_terms = enr.results[enr.results['Adjusted P-value'] < 0.05]
    if len(sig_terms) > 0:
        clusters_with_terms += 1
    print(f'Cluster {cl_id}: {len(sig_terms)} significant term(s)')

print(f'Clusters with significant terms: {clusters_with_terms} / {N_CLUSTERS}')
if clusters_with_terms < 3:
    print('WARNING: Few clusters enriched. On real data, check gene ID mapping or thresholds.')

shutil.rmtree(workdir, ignore_errors=True)   # remove all pipeline outputs; leave no stray files
print(f'\nPipeline complete: {len(temporal_genes)} temporal genes, {N_CLUSTERS} clusters, {clusters_with_terms} enriched')
