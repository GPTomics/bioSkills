import numpy as np
import matplotlib.pyplot as plt
from tslearn.clustering import TimeSeriesKMeans
from tslearn.preprocessing import TimeSeriesScalerMeanVariance
from sklearn.metrics import silhouette_score

np.random.seed(42)

# --- Simulate temporal expression profiles ---
# 8 timepoints: typical early-response time course
n_timepoints = 8
n_genes = 400
n_patterns = 4

patterns = np.array([
    [0, 1, 2, 2.5, 2, 1.5, 1, 0.5],     # early up
    [0, 0, 0.2, 0.5, 1, 2, 2.5, 3],      # late up
    [0, 0.5, 2, 3, 2, 0.5, 0, 0],        # transient
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

# --- Standardize for clustering ---
expr_3d = expr_mat[:, :, np.newaxis]
scaler = TimeSeriesScalerMeanVariance()
expr_scaled = scaler.fit_transform(expr_3d)

# --- Select k via silhouette score ---
# k range 2-10: below 2 is trivial; above 10 rarely adds meaning for 8 timepoints
k_range = range(2, 11)
sil_scores = []
for k in k_range:
    # metric='softdtw': differentiable approximation of DTW, faster for large datasets
    # gamma_sdtw=0.5: smoothing parameter; lower = closer to true DTW, higher = smoother
    model = TimeSeriesKMeans(n_clusters=k, metric='softdtw', gamma_sdtw=0.5, max_iter=30, random_state=42)
    labels = model.fit_predict(expr_scaled)
    # Euclidean silhouette as computational shortcut; DTW silhouette is O(n^2 * T^2)
    sil = silhouette_score(expr_scaled.squeeze(), labels, metric='euclidean')
    sil_scores.append(sil)

best_k = list(k_range)[np.argmax(sil_scores)]
print(f'Best k by silhouette: {best_k} (score: {max(sil_scores):.3f})')

# --- Cluster with optimal k ---
model = TimeSeriesKMeans(n_clusters=best_k, metric='softdtw', gamma_sdtw=0.5, max_iter=50, random_state=42)
labels = model.fit_predict(expr_scaled)

for k in range(best_k):
    print(f'  Cluster {k}: {np.sum(labels == k)} genes')

# --- Visualization ---
fig, axes = plt.subplots(2, 3, figsize=(15, 8))

axes[0, 0].plot(list(k_range), sil_scores, 'o-', color='steelblue', linewidth=2)
axes[0, 0].axvline(best_k, color='red', linestyle='--', alpha=0.7)
axes[0, 0].set_xlabel('Number of clusters (k)')
axes[0, 0].set_ylabel('Silhouette score')
axes[0, 0].set_title('Cluster selection')

timepoint_labels = ['0h', '2h', '4h', '8h', '12h', '24h', '36h', '48h']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
for k in range(min(best_k, 5)):
    ax = axes[(k + 1) // 3, (k + 1) % 3]
    cluster_mask = labels == k
    cluster_profiles = expr_scaled[cluster_mask].squeeze()

    for profile in cluster_profiles:
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
plt.savefig('tslearn_dtw_clusters.png', dpi=150, bbox_inches='tight')
plt.show()
print('\nPlot saved to tslearn_dtw_clusters.png')

# --- Export cluster assignments ---
cluster_df_rows = []
for i in range(n_genes):
    cluster_df_rows.append(f'gene_{i},{labels[i]}')

with open('tslearn_cluster_assignments.csv', 'w') as f:
    f.write('gene,cluster\n')
    f.write('\n'.join(cluster_df_rows))

print('Cluster assignments saved to tslearn_cluster_assignments.csv')
