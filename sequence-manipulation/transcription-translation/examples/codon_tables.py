'''Using alternative codon tables for translation'''
from Bio.Seq import Seq
from Bio.Data import CodonTable

# Standard translation
print('=== Standard Codon Table (1) ===')
seq = Seq('ATGTTTGGTAGA')
standard = seq.translate()
print(f'Sequence: {seq}')
print(f'Standard: {standard}')

# Bacterial/Archaeal (table 11)
print('\n=== Bacterial Codon Table (11) ===')
bacterial = seq.translate(table=11)
print(f'Bacterial: {bacterial}')

# Vertebrate Mitochondrial (table 2)
print('\n=== Vertebrate Mitochondrial (2) ===')
mito_seq = Seq('ATGAGACGA')
mito = mito_seq.translate(table=2)
standard = mito_seq.translate(table=1)
print(f'Sequence: {mito_seq}')
print(f'Standard: {standard}')
print(f'Mitochondrial: {mito}')

# Inspect codon table
print('\n=== Codon Table Details ===')
table = CodonTable.unambiguous_dna_by_id[1]
print(f'Name: {table.names}')
print(f'Start codons: {table.start_codons}')
print(f'Stop codons: {table.stop_codons}')

# CDS validation
print('\n=== CDS Translation ===')
valid_cds = Seq('ATGTTTGGTTAA')
try:
    protein = valid_cds.translate(cds=True)
    print(f'CDS: {valid_cds}')
    print(f'Protein: {protein}')
except Exception as e:
    print(f'Error: {e}')
