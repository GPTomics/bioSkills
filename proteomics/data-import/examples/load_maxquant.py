'''Load and filter MaxQuant proteinGroups.txt output'''
import pandas as pd
import numpy as np

protein_groups = pd.read_csv('proteinGroups.txt', sep='\t', low_memory=False)
print(f'Loaded {len(protein_groups)} protein groups')

# Filter contaminants, reverse hits, and proteins identified only by site
mask = (
    (protein_groups['Potential contaminant'] != '+') &
    (protein_groups['Reverse'] != '+') &
    (protein_groups['Only identified by site'] != '+')
)
filtered = protein_groups[mask].copy()
print(f'After filtering: {len(filtered)} protein groups')

# Extract LFQ intensity columns
lfq_cols = [c for c in filtered.columns if c.startswith('LFQ intensity')]
print(f'Found {len(lfq_cols)} samples')

# Create clean intensity matrix
intensity_matrix = filtered[['Protein IDs', 'Gene names', 'Protein names'] + lfq_cols].copy()
intensity_matrix.columns = ['protein_ids', 'gene_names', 'protein_names'] + [c.replace('LFQ intensity ', '') for c in lfq_cols]

# Replace 0 with NaN (MaxQuant uses 0 for missing)
sample_cols = intensity_matrix.columns[3:]
intensity_matrix[sample_cols] = intensity_matrix[sample_cols].replace(0, np.nan)

# Log2 transform
intensity_matrix[sample_cols] = np.log2(intensity_matrix[sample_cols])

# Missing value summary
missing_pct = 100 * intensity_matrix[sample_cols].isna().sum().sum() / intensity_matrix[sample_cols].size
print(f'Missing values: {missing_pct:.1f}%')

intensity_matrix.to_csv('intensity_matrix.csv', index=False)
