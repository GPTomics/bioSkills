'''Proteomics QC analysis workflow'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

intensity_matrix = pd.read_csv('normalized_imputed.csv', index_col=0)
sample_info = pd.read_csv('sample_info.csv', index_col=0)

print('=== Proteomics QC Report ===\n')
print(f'Samples: {intensity_matrix.shape[1]}')
print(f'Proteins: {intensity_matrix.shape[0]}')

# Missing values
missing_pct = 100 * intensity_matrix.isna().sum().sum() / intensity_matrix.size
print(f'\nOverall missing: {missing_pct:.1f}%')

missing_per_sample = 100 * intensity_matrix.isna().sum() / len(intensity_matrix)
print(f'Missing per sample: {missing_per_sample.min():.1f}% - {missing_per_sample.max():.1f}%')

# Correlation analysis
corr_matrix = intensity_matrix.corr()
upper_tri = corr_matrix.values[np.triu_indices_from(corr_matrix, k=1)]
print(f'\nMedian pairwise correlation: {np.median(upper_tri):.3f}')
print(f'Min correlation: {np.min(upper_tri):.3f}')

# Identify potential outliers (low correlation)
mean_corr = corr_matrix.mean()
outliers = mean_corr[mean_corr < mean_corr.mean() - 2 * mean_corr.std()]
if len(outliers) > 0:
    print(f'\nPotential outliers (low correlation): {list(outliers.index)}')

# PCA for batch effects
imputed = intensity_matrix.fillna(intensity_matrix.median())
scaled = StandardScaler().fit_transform(imputed.T)
pca = PCA(n_components=5)
pcs = pca.fit_transform(scaled)

print(f'\nPCA variance explained:')
for i, var in enumerate(pca.explained_variance_ratio_[:3]):
    print(f'  PC{i+1}: {100*var:.1f}%')

# Plots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Correlation heatmap
sns.heatmap(corr_matrix, cmap='RdBu_r', center=1, vmin=0.8, vmax=1, ax=axes[0, 0], annot=False)
axes[0, 0].set_title('Sample Correlation')

# PCA
pc_df = pd.DataFrame({'PC1': pcs[:, 0], 'PC2': pcs[:, 1]}, index=intensity_matrix.columns)
if 'condition' in sample_info.columns:
    pc_df['condition'] = sample_info['condition']
    for cond in pc_df['condition'].unique():
        mask = pc_df['condition'] == cond
        axes[0, 1].scatter(pc_df.loc[mask, 'PC1'], pc_df.loc[mask, 'PC2'], label=cond, s=50)
    axes[0, 1].legend()
else:
    axes[0, 1].scatter(pc_df['PC1'], pc_df['PC2'], s=50)
axes[0, 1].set_xlabel(f'PC1 ({100*pca.explained_variance_ratio_[0]:.1f}%)')
axes[0, 1].set_ylabel(f'PC2 ({100*pca.explained_variance_ratio_[1]:.1f}%)')
axes[0, 1].set_title('PCA')

# Intensity distributions
for col in intensity_matrix.columns[:10]:
    axes[1, 0].hist(intensity_matrix[col].dropna(), bins=50, alpha=0.5, label=col)
axes[1, 0].set_xlabel('Log2 Intensity')
axes[1, 0].set_ylabel('Frequency')
axes[1, 0].set_title('Intensity Distributions')

# Missing values per sample
axes[1, 1].bar(range(len(missing_per_sample)), missing_per_sample.values)
axes[1, 1].axhline(30, color='red', linestyle='--', label='30% threshold')
axes[1, 1].set_xlabel('Sample')
axes[1, 1].set_ylabel('Missing %')
axes[1, 1].set_title('Missing Values per Sample')

plt.tight_layout()
plt.savefig('qc_report.pdf')
print('\nSaved: qc_report.pdf')

# QC summary
print('\n=== QC Summary ===')
flags = []
if missing_pct > 30:
    flags.append('HIGH_MISSING')
if np.min(upper_tri) < 0.8:
    flags.append('LOW_CORRELATION')
if len(outliers) > 0:
    flags.append('OUTLIERS_DETECTED')

if flags:
    print(f'Flags: {", ".join(flags)}')
else:
    print('Status: PASS - No QC issues detected')
