'''Basic codon usage analysis'''
from Bio.Seq import Seq
from Bio.SeqUtils import GC123
from Bio.Data import CodonTable
from collections import Counter

def count_codons(seq):
    seq_str = str(seq).upper()
    codons = [seq_str[i:i+3] for i in range(0, len(seq_str) - 2, 3)]
    return Counter(codons)

def codon_frequencies(seq):
    counts = count_codons(seq)
    total = sum(counts.values())
    return {codon: count / total for codon, count in counts.items()}

# Example coding sequence
seq = Seq('ATGCGATCGATCGATCGATCGATCGATCGATCGTAA')
print(f'Sequence: {seq}')
print(f'Length: {len(seq)} bp ({len(seq)//3} codons)')

# Count codons
print('\n=== Codon Counts ===')
counts = count_codons(seq)
for codon, count in sorted(counts.items()):
    print(f'{codon}: {count}')

# Frequencies
print('\n=== Codon Frequencies ===')
freqs = codon_frequencies(seq)
for codon, freq in sorted(freqs.items(), key=lambda x: x[1], reverse=True):
    bar = '#' * int(freq * 50)
    print(f'{codon}: {freq:.3f} {bar}')

# GC at codon positions
print('\n=== GC at Codon Positions ===')
gc_total, gc_pos1, gc_pos2, gc_pos3 = GC123(seq)
print(f'Total GC:        {gc_total:.1f}%')
print(f'1st position:    {gc_pos1:.1f}%')
print(f'2nd position:    {gc_pos2:.1f}%')
print(f'3rd position:    {gc_pos3:.1f}% (wobble)')

# Amino acid usage
print('\n=== Amino Acid Usage ===')
table = CodonTable.unambiguous_dna_by_id[1]
aa_counts = Counter()
for codon, count in counts.items():
    if codon in table.stop_codons:
        aa_counts['*'] += count
    elif codon in table.forward_table:
        aa_counts[table.forward_table[codon]] += count

for aa, count in sorted(aa_counts.items(), key=lambda x: x[1], reverse=True):
    print(f'{aa}: {count}')
