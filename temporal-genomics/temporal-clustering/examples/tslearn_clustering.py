# Reference: tslearn 0.8+, scikit-learn 1.4+ | Verify API if version differs
import os
import shutil
import tempfile
import numpy as np
import matplotlib.pyplot as plt
from tslearn.clustering import TimeSeriesKMeans, silhouette_score
from tslearn.preprocessing import TimeSeriesScalerMeanVariance

np.random.seed(42)

# --- Simulate temporal expression profiles (already selected as temporally variable) ---
# 8 timepoints: typical early-response time course
n_timepoints = 8
n_genes = 400
n_patterns = 4

patterns = np.array([
    [0, 1, 2, 2.5, 2, 1.5, 1, 0.5],          # early up
    [0, 0, 0.2, 0.5, 1, 2, 2.5, 3],          # late up
    [0, 0.5, 2, 3, 2, 0.5, 0, 0],            # transient
    [0, -0.5, -1, -1.5, -2, -2.5, -2, -1.5]  # down
])

expr_mat = np.zeros((n_genes, n_timepoints))
true_labels = np.zeros(n_genes, dtype=int)
for i in range(n_genes):
    pattern_idx = i % n_patterns
    true_labels[i] = pattern_idx
    scale = np.random.uniform(0.8, 1.5)
    base = np.random.uniform(6, 10)
    # SD = 0.3: moderate noise for normalized log-expression
    expr_mat[i] = base + patterns[pattern_idx] * scale + np.random.normal(0, 0.3, n_timepoints)

# --- Standardize (mandatory: z-score per gene so shape, not magnitude, drives clustering) ---
expr_scaled = TimeSeriesScalerMeanVariance().fit_transform(expr_mat[:, :, np.newaxis])

# --- Distance metric choice -------------------------------------------------------------------
# Default to Euclidean-on-zscore (phase-SENSITIVE, cheap, no fabricated structure). Escalate to
# DTW ONLY for real, expected phase shifts, and ALWAYS constrain the warping window: tslearn's
# default global_constraint=None is the singularity-prone config that can invent co-regulation.
# These simulated profiles have NO phase shift, so Euclidean is the correct choice here; the DTW
# block is shown constrained for the phase-shift case.
metric = 'euclidean'
# sakoe_chiba_radius=2: warping-window half-width in timepoints; small for tightly sampled data.
dtw_params = {'global_constraint': 'sakoe_chiba', 'sakoe_chiba_radius': 2}

# --- Select k, scoring under the SAME geometry that forms the clusters ---
# Scoring DTW clusters with a Euclidean silhouette is geometrically inconsistent and can mis-rank k.
# tslearn.clustering.silhouette_score takes metric='dtw'/'softdtw' and precomputes the matching
# distances internally, so k is ranked in the same space the clusters were built in.
# k range 2-10: below 2 is trivial; above 10 rarely adds meaning for 8 timepoints.
k_range = range(2, 11)
sil_scores = []
for k in k_range:
    if metric == 'dtw':
        model = TimeSeriesKMeans(n_clusters=k, metric='dtw', metric_params=dtw_params, max_iter=30, random_state=42)
        labels = model.fit_predict(expr_scaled)
        sil = silhouette_score(expr_scaled, labels, metric='dtw', metric_params=dtw_params)
    else:
        model = TimeSeriesKMeans(n_clusters=k, metric='euclidean', max_iter=30, random_state=42)
        labels = model.fit_predict(expr_scaled)
        sil = silhouette_score(expr_scaled, labels, metric='euclidean')
    sil_scores.append(sil)

best_k = list(k_range)[int(np.argmax(sil_scores))]
print(f'Best k by silhouette ({metric}): {best_k} (score: {max(sil_scores):.3f})')

# --- Cluster with chosen k ---
if metric == 'dtw':
    model = TimeSeriesKMeans(n_clusters=best_k, metric='dtw', metric_params=dtw_params, max_iter=50, random_state=42)
else:
    model = TimeSeriesKMeans(n_clusters=best_k, metric='euclidean', max_iter=50, random_state=42)
labels = model.fit_predict(expr_scaled)

for k in range(best_k):
    print(f'  Cluster {k}: {np.sum(labels == k)} genes')

# --- Visualization ---
fig, axes = plt.subplots(2, 3, figsize=(15, 8))

axes[0, 0].plot(list(k_range), sil_scores, 'o-', color='steelblue', linewidth=2)
axes[0, 0].axvline(best_k, color='red', linestyle='--', alpha=0.7)
axes[0, 0].set_xlabel('Number of clusters (k)')
axes[0, 0].set_ylabel(f'Silhouette score ({metric})')
axes[0, 0].set_title('Cluster selection')

timepoint_labels = ['0h', '2h', '4h', '8h', '12h', '24h', '36h', '48h']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
for k in range(min(best_k, 5)):
    ax = axes[(k + 1) // 3, (k + 1) % 3]
    cluster_mask = labels == k
    cluster_profiles = expr_scaled[cluster_mask].squeeze()

    for profile in np.atleast_2d(cluster_profiles):
        ax.plot(range(n_timepoints), profile, alpha=0.1, color=colors[k])

    centroid = model.cluster_centers_[k].squeeze()
    ax.plot(range(n_timepoints), centroid, color='black', linewidth=2.5)
    ax.set_xticks(range(n_timepoints))
    ax.set_xticklabels(timepoint_labels, rotation=45, fontsize=8)
    ax.set_title(f'Cluster {k} (n={np.sum(cluster_mask)})')
    ax.set_ylabel('Standardized expression')

for idx in range(min(best_k, 5) + 1, 6):
    axes[idx // 3, idx % 3].set_visible(False)

plt.tight_layout()

# Write outputs to a temp dir and remove it so the example leaves no stray files
out_dir = tempfile.mkdtemp(prefix='tslearn_clusters_')
plt.savefig(os.path.join(out_dir, 'tslearn_clusters.png'), dpi=150, bbox_inches='tight')

with open(os.path.join(out_dir, 'tslearn_cluster_assignments.csv'), 'w') as f:
    f.write('gene,cluster\n')
    f.write('\n'.join(f'gene_{i},{labels[i]}' for i in range(n_genes)))

print(f'\nPlot + cluster assignments written to {out_dir} (removed on exit)')
shutil.rmtree(out_dir)
