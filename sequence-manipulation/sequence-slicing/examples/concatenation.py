'''Joining and concatenating sequences'''
from Bio.Seq import Seq

# Basic concatenation
print('=== Basic Concatenation ===')
seq1 = Seq('ATGC')
seq2 = Seq('GGGG')
seq3 = Seq('TTAA')
combined = seq1 + seq2 + seq3
print(f'{seq1} + {seq2} + {seq3} = {combined}')

# Concatenate with string
print('\n=== Concatenate with String ===')
seq = Seq('ATGC')
extended = seq + 'NNNN' + seq
print(f'Extended: {extended}')

# Join multiple sequences
print('\n=== Join with Linker ===')
seqs = [Seq('ATGC'), Seq('GGGG'), Seq('TTAA')]
linker = Seq('NNN')
joined = linker.join(seqs)
print(f'Joined with NNN: {joined}')

# Simulate exon splicing
print('\n=== Exon Splicing ===')
genomic = Seq('ATGC' + 'nnnnnnnn' + 'GGGG' + 'nnnnnn' + 'TTAA')
exon1, exon2, exon3 = Seq('ATGC'), Seq('GGGG'), Seq('TTAA')
mrna = exon1 + exon2 + exon3
print(f'Genomic: {genomic}')
print(f'mRNA (exons joined): {mrna}')

# Build sequence from parts
print('\n=== Build from Parts ===')
parts = ['ATG', 'CGA', 'TCG', 'ATC', 'TAA']
full = Seq(''.join(parts))
print(f'Built from codons: {full}')
