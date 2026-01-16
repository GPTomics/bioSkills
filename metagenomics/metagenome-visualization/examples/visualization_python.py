import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

abundance = pd.read_csv('merged_abundance.txt', sep='\t', index_col=0)
species = abundance[abundance.index.str.contains('s__') & ~abundance.index.str.contains('t__')]
species.index = species.index.str.split('|').str[-1].str.replace('s__', '')

top_species = species.sum(axis=1).nlargest(10).index
species_top = species.loc[top_species].copy()
species_top.loc['Other'] = species.drop(top_species).sum()

fig, ax = plt.subplots(figsize=(12, 6))
species_top.T.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')
ax.set_xlabel('Sample')
ax.set_ylabel('Relative Abundance (%)')
ax.set_title('Species Composition')
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
plt.tight_layout()
plt.savefig('stacked_bar.png', dpi=300, bbox_inches='tight')
plt.close()

top20 = species.sum(axis=1).nlargest(20).index
fig, ax = plt.subplots(figsize=(12, 10))
sns.heatmap(species.loc[top20], cmap='YlOrRd', ax=ax, cbar_kws={'label': 'Abundance (%)'})
ax.set_xlabel('Sample')
ax.set_ylabel('Species')
ax.set_title('Top 20 Species Abundance Heatmap')
plt.tight_layout()
plt.savefig('heatmap.png', dpi=300)
plt.close()

species_t = species.T
scaler = StandardScaler()
scaled = scaler.fit_transform(species_t)
pca = PCA(n_components=2)
pca_result = pca.fit_transform(scaled)

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(pca_result[:, 0], pca_result[:, 1], s=100)
for i, sample in enumerate(species_t.index):
    ax.annotate(sample, (pca_result[i, 0] + 0.1, pca_result[i, 1] + 0.1))
ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
ax.set_title('PCA of Sample Composition')
plt.tight_layout()
plt.savefig('pca.png', dpi=300)
plt.close()

print('Visualizations saved: stacked_bar.png, heatmap.png, pca.png')
