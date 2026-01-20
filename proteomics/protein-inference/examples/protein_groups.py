'''Create protein groups from peptide-spectrum matches'''
import pandas as pd
from collections import defaultdict

psms = pd.read_csv('psms_1pct_fdr.tsv', sep='\t')
print(f'Loaded {len(psms)} PSMs')

# Build peptide-to-protein mapping
peptide_proteins = defaultdict(set)
for _, row in psms.iterrows():
    proteins = row['protein'].split(';')
    peptide_proteins[row['peptide']].update(proteins)

print(f'Unique peptides: {len(peptide_proteins)}')

# Build protein-to-peptides mapping
protein_peptides = defaultdict(set)
for peptide, proteins in peptide_proteins.items():
    for protein in proteins:
        protein_peptides[protein].add(peptide)

print(f'Proteins before inference: {len(protein_peptides)}')

# Count unique peptides per protein
protein_unique = {}
for protein, peptides in protein_peptides.items():
    unique_count = sum(1 for p in peptides if len(peptide_proteins[p]) == 1)
    protein_unique[protein] = unique_count

# Create protein groups (proteins with identical peptide sets)
peptide_set_groups = defaultdict(list)
for protein, peptides in protein_peptides.items():
    key = frozenset(peptides)
    peptide_set_groups[key].append(protein)

groups = []
for peptides, proteins in peptide_set_groups.items():
    proteins_sorted = sorted(proteins, key=lambda p: (-protein_unique.get(p, 0), p))
    group = {
        'protein_group': ';'.join(proteins_sorted),
        'lead_protein': proteins_sorted[0],
        'n_proteins': len(proteins),
        'peptides': ';'.join(sorted(peptides)),
        'n_peptides': len(peptides),
        'n_unique_peptides': sum(1 for p in peptides if len(peptide_proteins[p]) == 1)
    }
    groups.append(group)

groups_df = pd.DataFrame(groups)
groups_df = groups_df.sort_values('n_unique_peptides', ascending=False)

# Filter: require at least 2 peptides, 1 unique
confident = groups_df[(groups_df['n_peptides'] >= 2) & (groups_df['n_unique_peptides'] >= 1)]
print(f'\nProtein groups (>=2 peptides, >=1 unique): {len(confident)}')

confident.to_csv('protein_groups.csv', index=False)
print(f'Saved to protein_groups.csv')
