import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

counts = pd.read_csv('count_matrix.csv', index_col=0)
metadata = pd.read_csv('sample_info.csv', index_col=0)

print('=== RNA-seq QC Report ===\n')
print(f'Samples: {counts.shape[1]}')
print(f'Genes: {counts.shape[0]}')

print('\nLibrary sizes:')
lib_sizes = counts.sum()
print(lib_sizes.describe())

print('\nGenes detected per sample:')
genes_detected = (counts > 0).sum()
print(genes_detected.describe())

# Standard filtering: gene must have >=10 counts in >=3 samples
# min_counts=10: Removes genes with only 1-2 reads that could be noise
# min_samples=3: Ensures gene is expressed in multiple replicates, not just one outlier
# Expected: ~12,000-18,000 genes pass filtering for mammalian samples
min_counts, min_samples = 10, 3
gene_filter = (counts >= min_counts).sum(axis=1) >= min_samples
counts_filt = counts[gene_filter]
print(f'\nGenes after filtering: {len(counts_filt)}')

cpm = counts_filt * 1e6 / counts_filt.sum()
log_cpm = np.log2(cpm + 1)

fig, ax = plt.subplots(figsize=(10, 6))
sample_cor = log_cpm.corr()
sns.heatmap(sample_cor, annot=True, cmap='RdBu_r', center=0.95, vmin=0.8, vmax=1.0, ax=ax)
plt.title('Sample Correlation')
plt.tight_layout()
plt.savefig('qc_sample_correlation.png')
plt.close()

pca = PCA(n_components=2)
pca_result = pca.fit_transform(log_cpm.T)

fig, ax = plt.subplots(figsize=(8, 6))
for condition in metadata['condition'].unique():
    mask = metadata['condition'] == condition
    ax.scatter(pca_result[mask, 0], pca_result[mask, 1], label=condition, s=100)
    for i, name in enumerate(counts.columns[mask]):
        ax.annotate(name, (pca_result[mask, 0][i], pca_result[mask, 1][i]))
ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
ax.legend()
ax.set_title('PCA of Samples')
plt.tight_layout()
plt.savefig('qc_pca.png')
plt.close()

fig, ax = plt.subplots(figsize=(10, 4))
ax.bar(range(len(lib_sizes)), lib_sizes / 1e6)
ax.set_xticks(range(len(lib_sizes)))
ax.set_xticklabels(lib_sizes.index, rotation=45, ha='right')
ax.set_ylabel('Million Reads')
ax.set_title('Library Sizes')
plt.tight_layout()
plt.savefig('qc_library_sizes.png')
plt.close()

print('\nQC plots saved: qc_sample_correlation.png, qc_pca.png, qc_library_sizes.png')
