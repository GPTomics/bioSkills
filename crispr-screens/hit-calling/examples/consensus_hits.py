'''Consensus hit calling from multiple methods'''
import pandas as pd
import numpy as np

# Load results from different methods
mageck = pd.read_csv('mageck.gene_summary.txt', sep='\t')
bagel = pd.read_csv('bagel_bf.txt', sep='\t')

# Standardize column names
mageck = mageck[['id', 'neg|score', 'neg|fdr']].rename(columns={'id': 'gene'})
bagel = bagel[['Gene', 'BF']].rename(columns={'Gene': 'gene'})

# Merge
merged = mageck.merge(bagel, on='gene', how='outer')

# Call hits per method
merged['mageck_hit'] = merged['neg|fdr'] < 0.1
merged['bagel_hit'] = merged['BF'] > 5

# Consensus
merged['n_methods'] = merged['mageck_hit'].fillna(False).astype(int) + merged['bagel_hit'].fillna(False).astype(int)

# Results
print('=== Consensus Hit Calling ===')
print(f'MAGeCK hits (FDR<0.1): {merged["mageck_hit"].sum()}')
print(f'BAGEL2 hits (BF>5): {merged["bagel_hit"].sum()}')
print(f'Consensus hits (both): {(merged["n_methods"] == 2).sum()}')

# High confidence hits
high_conf = merged[merged['n_methods'] == 2].sort_values('neg|score')
print('\nTop consensus hits:')
print(high_conf[['gene', 'neg|score', 'neg|fdr', 'BF']].head(20).to_string(index=False))

# Save
high_conf.to_csv('consensus_hits.csv', index=False)
