'''Create UpSet plots with upsetplot for set intersection visualization'''

import matplotlib.pyplot as plt
from upsetplot import from_contents, UpSet, plot
import pandas as pd
import numpy as np

# --- ALTERNATIVE: Use real gene sets ---
# For realistic examples, load DE results or GO terms:
#
# de_results = pd.read_csv('de_results.csv')
# gene_sets = {
#     'Upregulated': de_results[de_results['log2FC'] > 1]['gene'].tolist(),
#     'Downregulated': de_results[de_results['log2FC'] < -1]['gene'].tolist(),
#     'Significant': de_results[de_results['padj'] < 0.05]['gene'].tolist()
# }

np.random.seed(42)

all_genes = [f'Gene{i}' for i in range(1, 501)]

# Simulated gene sets from different analyses
# Realistic scenario: DE genes across conditions
gene_sets = {
    'Treatment_A': np.random.choice(all_genes, 150, replace=False).tolist(),
    'Treatment_B': np.random.choice(all_genes, 130, replace=False).tolist(),
    'Timepoint_Early': np.random.choice(all_genes, 100, replace=False).tolist(),
    'Timepoint_Late': np.random.choice(all_genes, 180, replace=False).tolist(),
    'Pathway_Response': np.random.choice(all_genes, 90, replace=False).tolist()
}

# Add core overlap genes for biological realism
# ~20-30 genes often appear across multiple related conditions
core_genes = np.random.choice(all_genes, 25, replace=False).tolist()
for key in ['Treatment_A', 'Treatment_B', 'Pathway_Response']:
    gene_sets[key] = list(set(gene_sets[key] + core_genes))

# Convert to upsetplot format
data = from_contents(gene_sets)

# Basic UpSet plot
# show_counts displays intersection sizes on bars
fig, ax = plt.subplots(figsize=(12, 8))
plot(data, show_counts=True, fig=fig)
plt.savefig('upset_basic.png', dpi=150, bbox_inches='tight')
plt.close()

# Customized UpSet plot
upset = UpSet(data,
              subset_size='count',
              show_counts=True,
              show_percentages=False,  # Can be True to show % of total
              sort_by='cardinality',   # 'cardinality' = frequency, 'degree' = num sets
              sort_categories_by='cardinality',
              facecolor='#4DBBD5',
              element_size=46,
              intersection_plot_elements=15)  # Max intersections to show

fig = plt.figure(figsize=(14, 8))
upset.plot(fig=fig)
plt.suptitle('Gene Set Intersections', fontsize=14, y=1.02)
plt.savefig('upset_customized.png', dpi=300, bbox_inches='tight')
plt.close()

# UpSet with metadata
# Add attributes to visualize per-intersection statistics
np.random.seed(42)
n_elements = len(data)
# Simulated log fold changes and p-values
log2fc = np.random.normal(0, 1.5, n_elements)
pvalues = 10 ** np.random.uniform(-5, -0.5, n_elements)

df_with_attrs = data.to_frame()
df_with_attrs['log2FC'] = log2fc
df_with_attrs['pvalue'] = pvalues
df_with_attrs['significant'] = df_with_attrs['pvalue'] < 0.05

# Recreate multi-index for upsetplot
df_indexed = df_with_attrs.set_index(list(gene_sets.keys()))

# UpSet with category plot (boxplot of log2FC per intersection)
upset = UpSet(df_indexed, subset_size='count', show_counts=True)
# add_catplot shows attribute distribution per intersection
# kind options: 'box', 'violin', 'strip', 'swarm'
upset.add_catplot(value='log2FC', kind='box', color='#E64B35')

fig = plt.figure(figsize=(14, 10))
upset.plot(fig=fig)
plt.savefig('upset_with_boxplot.png', dpi=300, bbox_inches='tight')
plt.close()

print('UpSet plots saved: upset_basic.png, upset_customized.png, upset_with_boxplot.png')

# Print statistics
print('\nSet sizes:')
for name, genes in gene_sets.items():
    print(f'  {name}: {len(genes)} genes')
