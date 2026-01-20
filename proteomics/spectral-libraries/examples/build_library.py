'''Build and validate a spectral library from search results'''
import pandas as pd
import matplotlib.pyplot as plt

# Load search results (MaxQuant evidence.txt format)
evidence = pd.read_csv('evidence.txt', sep='\t')

# Filter for high-quality PSMs
high_quality = evidence[
    (evidence['Score'] > 80) &
    (evidence['PEP'] < 0.01) &
    (evidence['Reverse'] != '+') &
    (evidence['Potential contaminant'] != '+')
].copy()

# Select best spectrum per precursor (highest score)
best_psms = high_quality.sort_values('Score', ascending=False)
best_psms = best_psms.drop_duplicates(subset=['Modified sequence', 'Charge'])

print(f'High-quality precursors: {len(best_psms)}')
print(f'Unique proteins: {best_psms["Proteins"].nunique()}')

# Create library format (simplified)
library = pd.DataFrame({
    'ModifiedSequence': best_psms['Modified sequence'],
    'PrecursorCharge': best_psms['Charge'],
    'PrecursorMz': best_psms['m/z'],
    'NormalizedRetentionTime': best_psms['Retention time'] / best_psms['Retention time'].max(),
    'ProteinId': best_psms['Proteins'],
    'GeneName': best_psms['Gene names']
})

# Save library
library.to_csv('empirical_library.tsv', sep='\t', index=False)

# QC plots
fig, axes = plt.subplots(1, 3, figsize=(12, 4))

# RT distribution
axes[0].hist(library['NormalizedRetentionTime'], bins=50)
axes[0].set_xlabel('Normalized RT')
axes[0].set_ylabel('Precursors')
axes[0].set_title('RT Distribution')

# Charge distribution
charge_counts = library['PrecursorCharge'].value_counts().sort_index()
axes[1].bar(charge_counts.index, charge_counts.values)
axes[1].set_xlabel('Charge State')
axes[1].set_ylabel('Precursors')
axes[1].set_title('Charge Distribution')

# m/z distribution
axes[2].hist(library['PrecursorMz'], bins=50)
axes[2].set_xlabel('Precursor m/z')
axes[2].set_ylabel('Precursors')
axes[2].set_title('m/z Distribution')

plt.tight_layout()
plt.savefig('library_qc.png', dpi=150)
print('Library QC plot saved to library_qc.png')
