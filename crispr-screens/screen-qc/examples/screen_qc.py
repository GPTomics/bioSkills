'''CRISPR screen quality control'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load counts (MAGeCK format)
counts = pd.read_csv('screen.count.txt', sep='\t', index_col=0)
genes = counts['Gene']
count_matrix = counts.drop('Gene', axis=1)

print('=== CRISPR Screen QC Report ===\n')
print(f'sgRNAs: {len(count_matrix)}')
print(f'Genes: {genes.nunique()}')
print(f'Samples: {list(count_matrix.columns)}\n')

# Zero counts - <1% dropout indicates good library coverage; >5% suggests sequencing depth issues
print('--- Library Representation ---')
zero_pct = (count_matrix == 0).sum() / len(count_matrix) * 100
for sample, pct in zero_pct.items():
    status = 'PASS' if pct < 1 else 'WARN' if pct < 5 else 'FAIL'
    print(f'{sample}: {pct:.2f}% zero counts [{status}]')

# Gini index
def gini_index(x):
    x = np.sort(x[x > 0])
    n = len(x)
    return (n + 1 - 2 * np.sum(np.cumsum(x)) / np.cumsum(x)[-1]) / n

# Gini index measures count inequality; <0.2 is uniform, >0.3 indicates bias/bottleneck
print('\n--- Read Distribution (Gini) ---')
for sample in count_matrix.columns:
    g = gini_index(count_matrix[sample].values)
    status = 'PASS' if g < 0.2 else 'WARN' if g < 0.3 else 'FAIL'
    print(f'{sample}: {g:.3f} [{status}]')

# Correlation
# Replicate correlation: r>0.8 expected for good replicates; r<0.6 indicates technical issues
print('\n--- Replicate Correlation ---')
log_counts = np.log10(count_matrix + 1)
corr = log_counts.corr()
for i, c1 in enumerate(count_matrix.columns):
    for c2 in count_matrix.columns[i+1:]:
        r = corr.loc[c1, c2]
        status = 'PASS' if r > 0.8 else 'WARN' if r > 0.6 else 'FAIL'
        print(f'{c1} vs {c2}: r={r:.3f} [{status}]')

# Plot
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for sample in count_matrix.columns:
    axes[0].hist(np.log10(count_matrix[sample] + 1), bins=50, alpha=0.5, label=sample)
axes[0].set_xlabel('Log10(counts + 1)')
axes[0].set_ylabel('sgRNAs')
axes[0].legend()

import seaborn as sns
sns.heatmap(corr, annot=True, cmap='RdYlBu_r', vmin=0.5, vmax=1, ax=axes[1])
axes[1].set_title('Sample Correlation')

plt.tight_layout()
plt.savefig('screen_qc.png', dpi=150)
print('\nQC plots saved to screen_qc.png')
