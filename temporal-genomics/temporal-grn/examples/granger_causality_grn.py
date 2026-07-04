# Reference: statsmodels 0.14+, numpy 1.26+, pandas 2.2+ | Verify API if version differs
import os
import io
import tempfile
import contextlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statsmodels.tsa.api import VAR
from statsmodels.tsa.stattools import grangercausalitytests, adfuller
from statsmodels.stats.multitest import multipletests

np.random.seed(42)

# --- Simulate TF-target time-series with a causal lag ---
# 20 timepoints: enough that maxlag up to ~3 stays above the n > 3*maxlag+1 DoF floor.
# Genome time courses usually have 6-12; 20 here so the demo has any power at all.
n_timepoints = 20
n_tfs = 5
n_targets = 15
n_genes = n_tfs + n_targets
gene_names = [f'TF{i+1}' for i in range(n_tfs)] + [f'target{i+1}' for i in range(n_targets)]
expr_mat = np.zeros((n_genes, n_timepoints))

# AR(1) stationary sd = sigma_e / sqrt(1 - phi^2); phi=0.7, sigma_e=0.5 seeds t0 at the
# stationary variance so the series does not have a burn-in transient.
ar_sd = 0.5 / np.sqrt(1 - 0.7**2)

for i in range(n_tfs):
    expr_mat[i, 0] = np.random.normal(0, ar_sd)
    for t in range(1, n_timepoints):
        # AR(1) coefficient 0.7: persistent but mean-reverting (stationary regulator)
        expr_mat[i, t] = 0.7 * expr_mat[i, t - 1] + np.random.normal(0, 0.5)

n_true_edges = 10
true_edges = []
for j in range(n_targets):
    target_idx = n_tfs + j
    expr_mat[target_idx, 0] = np.random.normal(0, ar_sd)
    if j < n_true_edges:
        causal_tf = j % n_tfs
        true_edges.append((gene_names[causal_tf], gene_names[target_idx]))
        for t in range(1, n_timepoints):
            # 0.4 cross-coefficient: moderate regulation, detectable above noise=0.3
            expr_mat[target_idx, t] = (0.5 * expr_mat[target_idx, t - 1] +
                                       0.4 * expr_mat[causal_tf, t - 1] +
                                       np.random.normal(0, 0.3))
    else:
        for t in range(1, n_timepoints):
            expr_mat[target_idx, t] = 0.7 * expr_mat[target_idx, t - 1] + np.random.normal(0, 0.5)

expr_df = pd.DataFrame(expr_mat, index=gene_names, columns=[f't{i}' for i in range(n_timepoints)])

# --- Stationarity: decide, do not reflexively difference ---
# ADF null = unit root (non-stationary); p < 0.05 => stationary. ADF is nearly powerless at
# 20 timepoints, so its result is ADVISORY here: it flags these stationary-by-construction
# AR(1) series as non-stationary purely from low power. Differencing is NOT free -- it removes
# the trend that carries the regulatory signal, and over-differencing a stationary series
# destroys the lag structure. So model on LEVELS (the DGP is stationary); genuinely
# non-stationary data would be differenced UNIFORMLY across all genes (never per-gene, which
# corrupts the VAR F-test reference distribution).
nonstationary_frac = np.mean([adfuller(expr_df.loc[g])[1] > 0.05 for g in gene_names])
print(f'Fraction of genes ADF-flagged non-stationary: {nonstationary_frac:.2f} '
      f'(advisory only -- ADF is underpowered at {n_timepoints} timepoints)')
expr_model = expr_df  # stationary AR(1) by construction; do not over-difference

# --- Pairwise Granger, BIC-selected lag, single test per pair ---
# maxlag=3 is the search ceiling; BIC picks ONE lag per pair so the reported p-value is not
# the optimistic best-of-several (a within-pair multiple test). NOTE: verbose is dropped --
# it is deprecated since statsmodels 0.14; the function still prints by default in 0.14, so
# the call is wrapped to silence that without passing the deprecated argument.
maxlag = 3
tf_names = [f'TF{i+1}' for i in range(n_tfs)]
target_names = [f'target{i+1}' for i in range(n_targets)]


def granger_pvalue(pair_data, maxlag):
    # column 0 = response Y (target), column 1 = predictor X (TF): tests X -> Y.
    lag = max(1, int(VAR(pair_data).select_order(maxlag).bic))
    with contextlib.redirect_stdout(io.StringIO()):
        res = grangercausalitytests(pair_data, maxlag=[lag])  # list -> only this lag
    return res[lag][0]['ssr_ftest'][1], res[lag][0]['ssr_ftest'][0], lag


records = []
for tf in tf_names:
    for target in target_names:
        pair = np.column_stack([expr_model.loc[target].values, expr_model.loc[tf].values])
        p, f_stat, lag = granger_pvalue(pair, maxlag)
        records.append({'tf': tf, 'target': target, 'p_value': p, 'f_stat': f_stat, 'lag': lag})

results_df = pd.DataFrame(records)

# --- Multiple testing across pairs ---
# multipletests default is Holm-Sidak, NOT BH; force fdr_bh explicitly.
results_df['q_value'] = multipletests(results_df['p_value'], method='fdr_bh')[1]

# q < 0.05: conventional FDR threshold for edge significance
significant = results_df[results_df['q_value'] < 0.05].sort_values('q_value')
print(f'Significant edges (q < 0.05): {len(significant)} / {len(results_df)}')
print(f'True edges in dataset: {len(true_edges)}')

# --- Evaluate against ground truth ---
predicted_edges = set(zip(significant['tf'], significant['target']))
true_edge_set = set(true_edges)
tp = len(predicted_edges & true_edge_set)
fp = len(predicted_edges - true_edge_set)
fn = len(true_edge_set - predicted_edges)
precision = tp / (tp + fp) if (tp + fp) > 0 else 0
recall = tp / (tp + fn) if (tp + fn) > 0 else 0
print(f'Precision: {precision:.2f}, Recall: {recall:.2f} '
      f'(low precision is expected -- edges are hypotheses, not validated regulation)')

# --- Build adjacency matrix ---
all_genes_sorted = tf_names + target_names
adj_matrix = pd.DataFrame(0.0, index=all_genes_sorted, columns=all_genes_sorted)
for _, row in significant.iterrows():
    # -log10(q): edge weight; higher = stronger evidence
    adj_matrix.loc[row['tf'], row['target']] = -np.log10(row['q_value'])

# --- Visualization (written to a temp dir so no stray files land in the repo) ---
out_path = os.path.join(tempfile.gettempdir(), 'granger_grn_results.png')
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

axes[0].hist(results_df['p_value'], bins=20, color='steelblue', edgecolor='black')
axes[0].axvline(0.05, color='red', linestyle='--', label='p = 0.05')
axes[0].set_xlabel('P-value')
axes[0].set_ylabel('Count')
axes[0].set_title('Granger p-value distribution')
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
    axes[1].set_title('Top edges (green = true)')

im = axes[2].imshow(adj_matrix.values[:n_tfs, n_tfs:], cmap='YlOrRd', aspect='auto')
axes[2].set_xticks(range(n_targets))
axes[2].set_xticklabels(target_names, rotation=90, fontsize=7)
axes[2].set_yticks(range(n_tfs))
axes[2].set_yticklabels(tf_names, fontsize=8)
axes[2].set_title('TF-target adjacency (-log10 q)')
plt.colorbar(im, ax=axes[2], shrink=0.8)

plt.tight_layout()
plt.savefig(out_path, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Network figure rendered to {out_path} (removed on exit)')
os.remove(out_path)
