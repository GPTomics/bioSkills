'''Create customized volcano plots with matplotlib'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text

# --- ALTERNATIVE: Use real DE results ---
# For realistic volcanos, load DESeq2/edgeR output:
#
# df = pd.read_csv('deseq2_results.csv')
# df = df.dropna(subset=['padj'])

np.random.seed(42)

# Simulate DE results with realistic distribution
n_genes = 5000

# Most genes no change, subset up/down regulated
log2fc = np.concatenate([
    np.random.normal(2, 0.5, 200),      # Upregulated
    np.random.normal(-1.8, 0.6, 150),   # Downregulated
    np.random.normal(0, 0.5, 4650)      # No change
])

# P-values correlate with effect size (realistic)
pvalue = 10 ** (-np.abs(log2fc) * np.random.uniform(0.5, 3, n_genes) -
                np.random.uniform(0, 2, n_genes))
pvalue = np.clip(pvalue, 0, 1)

df = pd.DataFrame({
    'gene': [f'Gene{i}' for i in range(1, n_genes + 1)],
    'log2FoldChange': log2fc,
    'pvalue': pvalue
})

# BH adjustment for multiple testing
# Using statsmodels for broader compatibility (scipy.stats.false_discovery_control requires scipy 1.11+)
from statsmodels.stats.multitest import multipletests
_, df['padj'], _, _ = multipletests(df['pvalue'], method='fdr_bh')

# Thresholds
# FC > 1: 2-fold change, standard for RNA-seq
# padj < 0.05: 5% FDR, standard threshold
fc_threshold = 1
padj_threshold = 0.05

# Significance categories
def get_significance(row):
    if row['padj'] < padj_threshold and row['log2FoldChange'] > fc_threshold:
        return 'Up'
    elif row['padj'] < padj_threshold and row['log2FoldChange'] < -fc_threshold:
        return 'Down'
    return 'NS'

df['significance'] = df.apply(get_significance, axis=1)

print('Significant genes:')
print(f"  Upregulated: {sum(df['significance'] == 'Up')}")
print(f"  Downregulated: {sum(df['significance'] == 'Down')}")
print(f"  Not significant: {sum(df['significance'] == 'NS')}")

# Color mapping
color_map = {'Up': '#E64B35', 'Down': '#4DBBD5', 'NS': 'lightgray'}
colors = df['significance'].map(color_map)

# Basic volcano plot
fig, ax = plt.subplots(figsize=(8, 6))

ax.scatter(df['log2FoldChange'], -np.log10(df['pvalue']),
           c=colors, alpha=0.6, s=15, edgecolors='none')

# Threshold lines
ax.axhline(-np.log10(0.05), color='gray', linestyle='--', linewidth=1)
ax.axvline(-fc_threshold, color='gray', linestyle='--', linewidth=1)
ax.axvline(fc_threshold, color='gray', linestyle='--', linewidth=1)

ax.set_xlabel(r'$log_2$ Fold Change', fontsize=12)
ax.set_ylabel(r'$-log_{10}$ p-value', fontsize=12)
ax.set_title('Differential Expression', fontsize=14)

# Clean up axes
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('volcano_basic.png', dpi=300, bbox_inches='tight')
plt.close()

# Volcano with gene labels
# Label top 15 significant genes by p-value
sig_df = df[df['significance'] != 'NS'].nsmallest(15, 'pvalue')

fig, ax = plt.subplots(figsize=(10, 8))

ax.scatter(df['log2FoldChange'], -np.log10(df['pvalue']),
           c=colors, alpha=0.5, s=15, edgecolors='none')

# Add labels with adjust_text to prevent overlaps
texts = []
for _, row in sig_df.iterrows():
    texts.append(ax.text(row['log2FoldChange'],
                         -np.log10(row['pvalue']),
                         row['gene'],
                         fontsize=8))

# adjust_text repositions labels to avoid overlaps
# arrowprops draws connecting lines
adjust_text(texts,
            arrowprops=dict(arrowstyle='-', color='gray', lw=0.5),
            expand_points=(1.5, 1.5))

ax.axhline(-np.log10(0.05), color='gray', linestyle='--', linewidth=1)
ax.axvline(-1, color='gray', linestyle='--', linewidth=1)
ax.axvline(1, color='gray', linestyle='--', linewidth=1)

ax.set_xlabel(r'$log_2$ Fold Change', fontsize=12)
ax.set_ylabel(r'$-log_{10}$ p-value', fontsize=12)
ax.set_title('DE Genes with Labels', fontsize=14)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('volcano_labeled.png', dpi=300, bbox_inches='tight')
plt.close()

# Volcano highlighting specific genes
genes_of_interest = ['Gene1', 'Gene201', 'Gene350', 'Gene25', 'Gene180']
highlight_df = df[df['gene'].isin(genes_of_interest)]

fig, ax = plt.subplots(figsize=(10, 8))

ax.scatter(df['log2FoldChange'], -np.log10(df['pvalue']),
           c=colors, alpha=0.3, s=15, edgecolors='none')

# Highlight with different marker
ax.scatter(highlight_df['log2FoldChange'], -np.log10(highlight_df['pvalue']),
           c='yellow', s=80, edgecolors='black', linewidths=1, zorder=5)

# Labels for highlighted genes
texts = []
for _, row in highlight_df.iterrows():
    texts.append(ax.text(row['log2FoldChange'],
                         -np.log10(row['pvalue']),
                         row['gene'],
                         fontsize=10,
                         fontweight='bold'))

adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', lw=1))

ax.set_xlabel(r'$log_2$ Fold Change', fontsize=12)
ax.set_ylabel(r'$-log_{10}$ p-value', fontsize=12)
ax.set_title('Genes of Interest Highlighted', fontsize=14)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('volcano_highlighted.png', dpi=300, bbox_inches='tight')
plt.close()

print('Volcano plots saved: volcano_basic.png, volcano_labeled.png, volcano_highlighted.png')
