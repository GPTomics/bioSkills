'''Advanced sequence property analysis: GC123, GC skew, amino acid codes'''
from Bio.Seq import Seq
from Bio.SeqUtils import GC123, GC_skew, nt_search, seq1, seq3

# GC at codon positions
print('=== GC at Codon Positions (GC123) ===')
coding_seq = Seq('ATGCGATCGATCGATCGATCGATCGATCGATCG')
print(f'Sequence: {coding_seq}')
gc_total, gc_pos1, gc_pos2, gc_pos3 = GC123(coding_seq)
print(f'Total GC:    {gc_total:.1f}%')
print(f'Position 1:  {gc_pos1:.1f}%')
print(f'Position 2:  {gc_pos2:.1f}%')
print(f'Position 3:  {gc_pos3:.1f}% (wobble)')

# GC skew analysis
print('\n=== GC Skew Analysis ===')
genome_seq = Seq('GCGCGCGCGCATATATATAT' * 5)
print(f'Sequence length: {len(genome_seq)} bp')
skew = GC_skew(genome_seq, window=20)
print(f'GC skew values (window=20): {[f"{s:.2f}" for s in skew[:5]]}...')
cumulative = sum(skew)
print(f'Cumulative skew: {cumulative:.2f}')

# IUPAC-aware search with nt_search
print('\n=== IUPAC Search (nt_search) ===')
seq = 'ATGCGATCGATCGATNGATCATGC'
pattern = 'GATNGATC'  # N matches any base
result = nt_search(seq, pattern)
print(f'Sequence: {seq}')
print(f'Pattern:  {pattern}')
print(f'Result: {result}')

# Amino acid code conversion
print('\n=== Amino Acid Code Conversion ===')
three_letter = 'MetAlaGlyTrpCysArg'
one_letter = seq1(three_letter)
print(f'3-letter: {three_letter}')
print(f'1-letter: {one_letter}')

back_to_three = seq3(one_letter)
with_separator = seq3(one_letter, join='-')
print(f'Back to 3-letter: {back_to_three}')
print(f'With separator: {with_separator}')

# Protein charge profile
print('\n=== Protein Charge at Different pH ===')
from Bio.SeqUtils.ProtParam import ProteinAnalysis
protein = ProteinAnalysis('MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIA')
print(f'Protein length: {len(protein.sequence)} aa')
for ph in [4.0, 5.0, 6.0, 7.0, 7.4, 8.0, 9.0, 10.0]:
    charge = protein.charge_at_pH(ph)
    bar = '+' * int(max(0, charge)) + '-' * int(max(0, -charge))
    print(f'pH {ph}: {charge:+.1f} {bar}')

# Flexibility profile (first 20 residues)
print('\n=== Flexibility Profile (first 20 residues) ===')
flex = protein.flexibility()
for i, f in enumerate(flex[:20]):
    bar = '#' * int(f * 10)
    print(f'{protein.sequence[i]}{i+1:3d}: {f:.3f} {bar}')
