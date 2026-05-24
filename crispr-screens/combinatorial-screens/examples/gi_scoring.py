# Reference: pandas 2.2+, numpy 1.26+, scipy 1.12+ | Verify API if version differs
#
# Genetic-interaction (GI) scoring for combinatorial CRISPR screens.
# Identifies synthetic-lethal and synthetic-rescue paralog pairs.

import pandas as pd
import numpy as np
from scipy.stats import zscore

# === INPUTS ===
# paired_lfc.tsv: cassette_id, gene_A, gene_B, lfc (paired double-KO LFC vs control)
# single_lfc.tsv: gene, lfc (single-KO LFC vs control)
paired_df = pd.read_csv('paired_lfc.tsv', sep='\t')
single_df = pd.read_csv('single_lfc.tsv', sep='\t')

# === STEP 1: BUILD SINGLE-GENE LOOKUP ===
single_lookup = dict(zip(single_df['gene'], single_df['lfc']))

# === STEP 2: AGGREGATE CASSETTE-LEVEL TO PAIR-LEVEL ===
# Inzolia has 4 cassettes per pair; aggregate mean LFC across replicates
pair_lfc = paired_df.groupby(['gene_A', 'gene_B']).agg(
    paired_lfc_mean=('lfc', 'mean'),
    paired_lfc_std=('lfc', 'std'),
    n_cassettes=('lfc', 'count')
).reset_index()

# === STEP 3: COMPUTE EXPECTED ADDITIVE FROM SINGLETONS ===
pair_lfc['single_A_lfc'] = pair_lfc['gene_A'].map(single_lookup)
pair_lfc['single_B_lfc'] = pair_lfc['gene_B'].map(single_lookup)
pair_lfc['expected_additive'] = pair_lfc['single_A_lfc'] + pair_lfc['single_B_lfc']

# Drop pairs where singleton LFC is missing
pair_lfc = pair_lfc.dropna(subset=['single_A_lfc', 'single_B_lfc'])

# === STEP 4: COMPUTE GI SCORE ===
# GI = observed_paired - expected_additive
# Negative GI = synthetic lethal (paired more depleted than expected)
# Positive GI = synthetic rescue (paired less depleted than expected)
pair_lfc['gi_score'] = pair_lfc['paired_lfc_mean'] - pair_lfc['expected_additive']

# === STEP 5: Z-NORMALIZE GI ===
pair_lfc['gi_z'] = zscore(pair_lfc['gi_score'])

# === STEP 6: CLASSIFY ===
pair_lfc['gi_class'] = np.where(pair_lfc['gi_z'] < -2, 'synthetic_lethal',
                                 np.where(pair_lfc['gi_z'] > 2, 'synthetic_rescue',
                                          'no_interaction'))

# === STEP 7: OUTPUT ===
# Synthetic-lethal candidates (drug-target combinations)
sl_pairs = pair_lfc[pair_lfc['gi_class'] == 'synthetic_lethal'].sort_values('gi_z')
print(f'Synthetic-lethal pairs (GI z <-2): {len(sl_pairs)}')
print(sl_pairs[['gene_A', 'gene_B', 'paired_lfc_mean', 'expected_additive',
                 'gi_score', 'gi_z']].head(20).to_string(index=False))

# Synthetic-rescue (compensatory pathways)
sr_pairs = pair_lfc[pair_lfc['gi_class'] == 'synthetic_rescue'].sort_values('gi_z', ascending=False)
print(f'\nSynthetic-rescue pairs (GI z >2): {len(sr_pairs)}')

# === EXPORT ===
pair_lfc.to_csv('gi_scores.tsv', sep='\t', index=False)
sl_pairs.to_csv('synthetic_lethal_pairs.tsv', sep='\t', index=False)
sr_pairs.to_csv('synthetic_rescue_pairs.tsv', sep='\t', index=False)

# === VALIDATION ===
# For top synthetic-lethal hits, cross-validate against:
# 1. Known paralog pairs (e.g., MAPK1/MAPK3, AKT1/AKT2)
# 2. Multiple cell lines (Dede 2020 reported 79% reproducibility across 3 lines)
# 3. Orthogonal modality (CRISPRi if originally Cas9, or vice versa)
# 4. Arrayed validation in target cell line
known_paralog_pairs = [('MAPK1', 'MAPK3'), ('AKT1', 'AKT2'),
                        ('PIK3CA', 'PIK3CB'), ('HSP90AA1', 'HSP90AB1')]
recovered = [(a, b) for a, b in known_paralog_pairs
              if ((pair_lfc['gene_A'] == a) & (pair_lfc['gene_B'] == b) &
                  (pair_lfc['gi_class'] == 'synthetic_lethal')).any()
              or ((pair_lfc['gene_A'] == b) & (pair_lfc['gene_B'] == a) &
                  (pair_lfc['gi_class'] == 'synthetic_lethal')).any()]
print(f'\nKnown paralog pairs recovered: {len(recovered)}/{len(known_paralog_pairs)}')
