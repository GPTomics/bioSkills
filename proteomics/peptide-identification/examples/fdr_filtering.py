'''Calculate FDR and filter PSMs using target-decoy approach'''
import pandas as pd
import numpy as np

psms = pd.read_csv('search_results.tsv', sep='\t')
print(f'Loaded {len(psms)} PSMs')

# Identify decoy hits (common prefixes: REV_, DECOY_, XXX_)
psms['is_decoy'] = psms['protein'].str.startswith(('REV_', 'DECOY_', 'XXX_'))
n_targets = (~psms['is_decoy']).sum()
n_decoys = psms['is_decoy'].sum()
print(f'Targets: {n_targets}, Decoys: {n_decoys}')

# Sort by score (descending - higher is better)
psms_sorted = psms.sort_values('score', ascending=False).reset_index(drop=True)

# Calculate cumulative FDR
target_count = (~psms_sorted['is_decoy']).cumsum()
decoy_count = psms_sorted['is_decoy'].cumsum()
psms_sorted['fdr'] = decoy_count / target_count

# Calculate q-value (minimum FDR at this score or better)
psms_sorted['qvalue'] = psms_sorted['fdr'][::-1].cummin()[::-1]

# Filter to 1% FDR
psms_filtered = psms_sorted[psms_sorted['qvalue'] <= 0.01].copy()
psms_filtered = psms_filtered[~psms_filtered['is_decoy']]
print(f'PSMs at 1% FDR: {len(psms_filtered)}')

# Unique peptides
unique_peptides = psms_filtered['peptide'].nunique()
print(f'Unique peptides: {unique_peptides}')

psms_filtered.to_csv('psms_1pct_fdr.tsv', sep='\t', index=False)
