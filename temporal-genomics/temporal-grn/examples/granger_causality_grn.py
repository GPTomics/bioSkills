import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import grangercausalitytests, adfuller
from statsmodels.stats.multitest import multipletests

np.random.seed(42)

# --- Simulate TF-target time-series with causal lag ---
# 20 timepoints: sufficient for Granger with maxlag=3 (need n > 3*maxlag)
n_timepoints = 20
n_tfs = 5
n_targets = 15
n_genes = n_tfs + n_targets

gene_names = [f'TF{i+1}' for i in range(n_tfs)] + [f'target{i+1}' for i in range(n_targets)]
expr_mat = np.zeros((n_genes, n_timepoints))

for i in range(n_tfs):
    expr_mat[i, 0] = np.random.uniform(6, 10)
    for t in range(1, n_timepoints):
        # AR(1) process: autoregressive with coefficient 0.7; persistent but mean-reverting
        expr_mat[i, t] = 0.7 * expr_mat[i, t - 1] + np.random.normal(0, 0.5)

n_true_edges = 10
true_edges = []
for j in range(n_targets):
    target_idx = n_tfs + j
    expr_mat[target_idx, 0] = np.random.uniform(6, 10)
    if j < n_true_edges:
        # Assign a causal TF with lag-1 influence
        causal_tf = j % n_tfs
        true_edges.append((gene_names[causal_tf], gene_names[target_idx]))
        for t in range(1, n_timepoints):
            # 0.4 coefficient: moderate regulatory influence (detectable above noise)
            expr_mat[target_idx, t] = (0.5 * expr_mat[target_idx, t - 1] +
                                        0.4 * expr_mat[causal_tf, t - 1] +
                                        np.random.normal(0, 0.3))
    else:
        for t in range(1, n_timepoints):
            expr_mat[target_idx, t] = 0.7 * expr_mat[target_idx, t - 1] + np.random.normal(0, 0.5)

expr_df = pd.DataFrame(expr_mat, index=gene_names, columns=[f't{i}' for i in range(n_timepoints)])

# --- Stationarity: apply uniform first differencing ---
# Uniform differencing avoids mixing differenced and non-differenced genes
# which would violate VAR model assumptions in Granger causality
expr_df = expr_df.diff(axis=1).iloc[:, 1:]

# --- Pairwise Granger causality tests ---
# maxlag=2: tests lag 1 and 2; lag 1 = one sampling interval delay
# Need n > 3 * maxlag timepoints for reliable estimation
maxlag = 2
tf_names = [f'TF{i+1}' for i in range(n_tfs)]
target_names = [f'target{i+1}' for i in range(n_targets)]

results = []
for tf in tf_names:
    for target in target_names:
        pair_data = np.column_stack([expr_df.loc[target].values, expr_df.loc[tf].values])
        gc_result = grangercausalitytests(pair_data, maxlag=maxlag, verbose=False)
        min_p = min(gc_result[lag][0]['ssr_ftest'][1] for lag in range(1, maxlag + 1))
        best_lag = min(range(1, maxlag + 1), key=lambda l: gc_result[l][0]['ssr_ftest'][1])
        f_stat = gc_result[best_lag][0]['ssr_ftest'][0]
        results.append({'tf': tf, 'target': target, 'p_value': min_p,
                         'f_stat': f_stat, 'best_lag': best_lag})

results_df = pd.DataFrame(results)

# --- Multiple testing correction ---
# BH FDR: standard for genome-wide regulatory edge testing
reject, qvals, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
results_df['q_value'] = qvals

# q < 0.05: standard FDR threshold
significant = results_df[results_df['q_value'] < 0.05].sort_values('q_value')
print(f'\nSignificant edges (q < 0.05): {len(significant)} / {len(results_df)}')
print(f'True edges in dataset: {len(true_edges)}')

# --- Evaluate against ground truth ---
predicted_edges = set(zip(significant['tf'], significant['target']))
true_edge_set = set(true_edges)
tp = len(predicted_edges & true_edge_set)
fp = len(predicted_edges - true_edge_set)
fn = len(true_edge_set - predicted_edges)
precision = tp / (tp + fp) if (tp + fp) > 0 else 0
recall = tp / (tp + fn) if (tp + fn) > 0 else 0
print(f'Precision: {precision:.2f}, Recall: {recall:.2f}')

# --- Build adjacency matrix ---
all_genes_sorted = tf_names + target_names
adj_matrix = pd.DataFrame(0.0, index=all_genes_sorted, columns=all_genes_sorted)
for _, row in significant.iterrows():
    # -log10(q): edge weight; higher = stronger evidence of regulation
    adj_matrix.loc[row['tf'], row['target']] = -np.log10(row['q_value'])

# --- Visualization ---
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

axes[0].hist(results_df['p_value'], bins=20, color='steelblue', edgecolor='black')
axes[0].axvline(0.05, color='red', linestyle='--', label='p = 0.05')
axes[0].set_xlabel('P-value')
axes[0].set_ylabel('Count')
axes[0].set_title('Granger causality p-value distribution')
axes[0].legend()

if len(significant) > 0:
    top_edges = significant.head(15)
    edge_labels = [f'{r["tf"]}->{r["target"]}' for _, r in top_edges.iterrows()]
    colors = ['green' if (r['tf'], r['target']) in true_edge_set else 'gray'
              for _, r in top_edges.iterrows()]
    axes[1].barh(range(len(edge_labels)), -np.log10(top_edges['q_value']),
                 color=colors, edgecolor='black')
    axes[1].set_yticks(range(len(edge_labels)))
    axes[1].set_yticklabels(edge_labels, fontsize=8)
    axes[1].set_xlabel('-log10(q-value)')
    axes[1].set_title('Top regulatory edges (green = true)')

im = axes[2].imshow(adj_matrix.values[:n_tfs, n_tfs:], cmap='YlOrRd', aspect='auto')
axes[2].set_xticks(range(n_targets))
axes[2].set_xticklabels(target_names, rotation=90, fontsize=7)
axes[2].set_yticks(range(n_tfs))
axes[2].set_yticklabels(tf_names, fontsize=8)
axes[2].set_title('TF-target adjacency (-log10 q)')
plt.colorbar(im, ax=axes[2], shrink=0.8)

plt.tight_layout()
plt.savefig('granger_grn_results.png', dpi=150, bbox_inches='tight')
plt.show()
print('\nPlot saved to granger_grn_results.png')
