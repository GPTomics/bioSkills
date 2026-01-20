'''Label-free quantification normalization'''
import pandas as pd
import numpy as np
from scipy import stats

intensity_matrix = pd.read_csv('intensity_matrix.csv', index_col=0)
sample_cols = intensity_matrix.columns.tolist()

# Already log2 transformed, handle missing values
print(f'Input shape: {intensity_matrix.shape}')
print(f'Missing values: {100 * intensity_matrix.isna().sum().sum() / intensity_matrix.size:.1f}%')

# Median centering normalization
sample_medians = intensity_matrix[sample_cols].median()
global_median = sample_medians.median()
normalized = intensity_matrix[sample_cols] - sample_medians + global_median
print(f'Sample medians before: {sample_medians.values}')
print(f'Sample medians after: {normalized.median().values}')

# Filter proteins with too many missing values (>50% per group)
valid_per_protein = normalized.notna().sum(axis=1)
min_valid = len(sample_cols) // 2
filtered = normalized[valid_per_protein >= min_valid]
print(f'Proteins after filtering: {len(filtered)}')

# Imputation: MinProb (left-censored Gaussian)
def minprob_impute(col, downshift=1.8, width=0.3):
    valid = col.dropna()
    if len(valid) == 0:
        return col
    mean_shift = valid.mean() - downshift * valid.std()
    std_shift = valid.std() * width
    imputed = col.copy()
    missing_mask = col.isna()
    imputed[missing_mask] = np.random.normal(mean_shift, std_shift, missing_mask.sum())
    return imputed

imputed = filtered.apply(minprob_impute)
print(f'Missing after imputation: {imputed.isna().sum().sum()}')

# Quality check: correlation heatmap values
corr_matrix = imputed.corr()
print(f'Mean pairwise correlation: {corr_matrix.values[np.triu_indices_from(corr_matrix, k=1)].mean():.3f}')

imputed.to_csv('normalized_imputed.csv')
