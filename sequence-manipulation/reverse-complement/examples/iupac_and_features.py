'''IUPAC ambiguity reverse complement, U/gap handling, and minus-strand feature extraction'''
# Reference: biopython 1.83+ | Verify API if version differs
from Bio.Seq import Seq, MutableSeq
from Bio.SeqFeature import SeqFeature, SimpleLocation

print('=== IUPAC ambiguity codes (Biopython handles all 15 + X) ===')
ambig = Seq('ACGTRYSWKMBDHVN')
print(f'Original: {ambig}')
print(f'RevComp:  {ambig.reverse_complement()}')
print(f'Self-complementary check (S,W,N): {Seq("SWN").complement()}')

print('\n=== Case-insensitive ===')
print(f"reverse_complement of 'atRY': {Seq('atRY').reverse_complement()}")

print('\n=== DNA mode maps U->A and emits T; RNA method keeps U ===')
print(f"Seq('ACGU').reverse_complement():     {Seq('ACGU').reverse_complement()}")
print(f"Seq('ACGU').reverse_complement_rna(): {Seq('ACGU').reverse_complement_rna()}")

print('\n=== Gaps and non-table chars pass through unchanged ===')
print(f"reverse_complement of 'ATG-CGA--TY': {Seq('ATG-CGA--TY').reverse_complement()}")

print('\n=== Protein trap: silent garbage, no warning ===')
protein = Seq('MAIVMGR')
print(f'Protein {protein} -> {protein.reverse_complement()} (meaningless; guard on molecule_type)')

print('\n=== In-place only works on MutableSeq ===')
m = MutableSeq('ATGC')
m.reverse_complement(inplace=True)
print(f'MutableSeq after inplace reverse_complement: {m}')

print('\n=== Minus-strand feature: extract() already reverse-complements ===')
parent = Seq('AAATGGGCCCTTTAAA')
feature = SeqFeature(SimpleLocation(3, 12, strand=-1), type='CDS')
cds = feature.extract(parent)
print(f'Parent:               {parent}')
print(f'extract() coding seq: {cds}  (correct, already RC for minus strand)')
print(f'Double-RC bug would give wrong strand: {cds.reverse_complement()}  (do NOT do this)')
