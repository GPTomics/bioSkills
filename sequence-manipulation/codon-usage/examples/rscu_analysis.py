'''RSCU (Relative Synonymous Codon Usage) analysis'''
from Bio.Seq import Seq
from Bio.Data import CodonTable
from collections import Counter

def count_codons(seq):
    seq_str = str(seq).upper()
    codons = [seq_str[i:i+3] for i in range(0, len(seq_str) - 2, 3)]
    return Counter(codons)

def calculate_rscu(seq, table_id=1):
    table = CodonTable.unambiguous_dna_by_id[table_id]
    counts = count_codons(seq)

    # Build back table: amino acid -> list of codons
    back_table = {}
    for codon, aa in table.forward_table.items():
        back_table.setdefault(aa, []).append(codon)

    rscu = {}
    for aa, codons in back_table.items():
        total = sum(counts.get(c, 0) for c in codons)
        n_synonymous = len(codons)
        expected = total / n_synonymous if n_synonymous > 0 else 0
        for codon in codons:
            observed = counts.get(codon, 0)
            rscu[codon] = observed / expected if expected > 0 else 0
    return rscu

# Example sequence
seq = Seq('ATGCTTCTGCTACTGCTTCTACTGCTGCTACTGTAA')
print(f'Sequence: {seq}')

# Calculate RSCU
print('\n=== RSCU Analysis ===')
rscu = calculate_rscu(seq)

# Get codon table for amino acid names
table = CodonTable.unambiguous_dna_by_id[1]

# Group by amino acid
back_table = {}
for codon, aa in table.forward_table.items():
    back_table.setdefault(aa, []).append(codon)

print('RSCU values (1.0 = no bias, >1 = preferred, <1 = avoided):')
print()

for aa in sorted(back_table.keys()):
    codons = back_table[aa]
    if len(codons) == 1:
        continue  # Skip amino acids with only one codon

    codon_info = [(c, rscu.get(c, 0)) for c in codons]
    non_zero = [c for c, r in codon_info if r > 0]
    if not non_zero:
        continue

    print(f'{aa}:')
    for codon, r in sorted(codon_info, key=lambda x: x[1], reverse=True):
        if r > 0:
            bar = '+' * int(r * 5) if r > 1 else '-' * int((1 - r) * 5)
            status = 'preferred' if r > 1.2 else 'avoided' if r < 0.8 else 'neutral'
            print(f'  {codon}: {r:.2f} {bar} ({status})')
    print()

# Find rare codons
print('=== Rare Codons (RSCU < 0.5) ===')
rare = {c: r for c, r in rscu.items() if 0 < r < 0.5}
counts = count_codons(seq)
for codon, r in sorted(rare.items(), key=lambda x: x[1]):
    aa = table.forward_table.get(codon, '?')
    print(f'{codon} ({aa}): RSCU={r:.2f}, count={counts.get(codon, 0)}')
