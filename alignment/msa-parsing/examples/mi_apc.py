'''Mutual information with average-product correction (MI-APC) for coevolution detection.

Dunn SD, Wahl LM, Gloor GB. 2008. Mutual information without the influence of phylogeny
or entropy dramatically improves residue contact prediction. Bioinformatics 24:333-340.
APC removes background signal from per-column entropy and phylogenetic structure.
The APC term excludes the i==j diagonal so the per-column mean reflects only
off-diagonal MI contributions, matching the Dunn et al definition.
'''
# Reference: biopython 1.83+, numpy 1.26+ | Verify API if version differs

from Bio import AlignIO
import numpy as np

def joint_counts(col_i, col_j):
    pairs = np.array([col_i, col_j])
    unique_i = np.unique(col_i)
    unique_j = np.unique(col_j)
    table = np.zeros((len(unique_i), len(unique_j)))
    idx_i = {r: k for k, r in enumerate(unique_i)}
    idx_j = {r: k for k, r in enumerate(unique_j)}
    for a, b in pairs.T:
        table[idx_i[a], idx_j[b]] += 1
    return table

def mi_matrix_apc(alignment, min_pairs=20):
    seq_array = np.array([list(str(r.seq)) for r in alignment])
    n_cols = seq_array.shape[1]
    mi = np.zeros((n_cols, n_cols))
    valid = np.zeros((n_cols, n_cols), dtype=bool)
    for i in range(n_cols):
        for j in range(i + 1, n_cols):
            col_i = seq_array[:, i]
            col_j = seq_array[:, j]
            mask = (col_i != '-') & (col_j != '-')
            if mask.sum() < min_pairs:
                continue
            table = joint_counts(col_i[mask], col_j[mask])
            joint_p = table / table.sum()
            margin_i = joint_p.sum(axis=1, keepdims=True)
            margin_j = joint_p.sum(axis=0, keepdims=True)
            with np.errstate(divide='ignore', invalid='ignore'):
                pmi = np.where(joint_p > 0, joint_p * np.log2(joint_p / (margin_i * margin_j)), 0)
            mi[i, j] = mi[j, i] = pmi.sum()
            valid[i, j] = valid[j, i] = True
    off_diag = ~np.eye(n_cols, dtype=bool) & valid
    column_means = np.zeros(n_cols)
    for k in range(n_cols):
        row_mask = off_diag[k]
        if row_mask.any():
            column_means[k] = mi[k, row_mask].mean()
    overall_mean = mi[off_diag].mean() if off_diag.any() else 0.0
    apc = np.outer(column_means, column_means) / overall_mean if overall_mean > 0 else np.zeros_like(mi)
    return mi - apc

if __name__ == '__main__':
    alignment = AlignIO.read('alignment.fasta', 'fasta')
    corrected = mi_matrix_apc(alignment)

    upper = np.triu_indices(corrected.shape[0], k=1)
    top_pairs = sorted(zip(corrected[upper], upper[0], upper[1]), reverse=True)[:20]
    print('Top 20 coevolving column pairs (MI-APC, bits):')
    for score, i, j in top_pairs:
        print(f'  {i:4d}-{j:4d}: {score:6.3f}')
