'''Codon optimization for heterologous expression'''
from Bio.Seq import Seq
from Bio.Data import CodonTable

# Preferred codons for common expression hosts
ECOLI_PREFERRED = {
    'A': 'GCG', 'R': 'CGT', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC',
    'Q': 'CAG', 'E': 'GAA', 'G': 'GGT', 'H': 'CAC', 'I': 'ATT',
    'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTC', 'P': 'CCG',
    'S': 'TCT', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAC', 'V': 'GTT',
    '*': 'TAA'
}

YEAST_PREFERRED = {
    'A': 'GCT', 'R': 'AGA', 'N': 'AAC', 'D': 'GAC', 'C': 'TGT',
    'Q': 'CAA', 'E': 'GAA', 'G': 'GGT', 'H': 'CAC', 'I': 'ATT',
    'L': 'TTG', 'K': 'AAG', 'M': 'ATG', 'F': 'TTC', 'P': 'CCA',
    'S': 'TCT', 'T': 'ACT', 'W': 'TGG', 'Y': 'TAC', 'V': 'GTT',
    '*': 'TAA'
}

def optimize_sequence(dna_seq, preferred_codons, table_id=1):
    '''Optimize codon usage while preserving amino acid sequence'''
    table = CodonTable.unambiguous_dna_by_id[table_id]
    seq_str = str(dna_seq).upper()

    optimized_codons = []
    changes = []

    for i in range(0, len(seq_str) - 2, 3):
        original_codon = seq_str[i:i+3]

        # Get amino acid
        if original_codon in table.stop_codons:
            aa = '*'
        elif original_codon in table.forward_table:
            aa = table.forward_table[original_codon]
        else:
            optimized_codons.append(original_codon)
            continue

        # Get preferred codon
        preferred = preferred_codons.get(aa, original_codon)
        optimized_codons.append(preferred)

        if preferred != original_codon:
            changes.append((i // 3, original_codon, preferred, aa))

    return Seq(''.join(optimized_codons)), changes

# Example: optimize a sequence for E. coli
original = Seq('ATGCTTCTGCTAAGCCGCTGA')
print('=== Original Sequence ===')
print(f'DNA: {original}')
protein = original.translate()
print(f'Protein: {protein}')

# Optimize for E. coli
print('\n=== E. coli Optimization ===')
ecoli_opt, ecoli_changes = optimize_sequence(original, ECOLI_PREFERRED)
print(f'Optimized: {ecoli_opt}')
print(f'Verifies: {ecoli_opt.translate()}')
print(f'Changes made: {len(ecoli_changes)}')
for pos, old, new, aa in ecoli_changes:
    print(f'  Position {pos}: {old} -> {new} ({aa})')

# Optimize for yeast
print('\n=== Yeast Optimization ===')
yeast_opt, yeast_changes = optimize_sequence(original, YEAST_PREFERRED)
print(f'Optimized: {yeast_opt}')
print(f'Verifies: {yeast_opt.translate()}')
print(f'Changes made: {len(yeast_changes)}')
for pos, old, new, aa in yeast_changes:
    print(f'  Position {pos}: {old} -> {new} ({aa})')

# Compare sequences
print('\n=== Comparison ===')
print(f'Original: {original}')
print(f'E. coli:  {ecoli_opt}')
print(f'Yeast:    {yeast_opt}')

# Show differences
print('\n=== Nucleotide Differences ===')
for i, (o, e, y) in enumerate(zip(str(original), str(ecoli_opt), str(yeast_opt))):
    if o != e or o != y:
        print(f'Position {i}: Original={o}, E.coli={e}, Yeast={y}')
