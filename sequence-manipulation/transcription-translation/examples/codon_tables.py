'''Standard vs mitochondrial translation and cds=True validation'''
# Reference: biopython 1.83+ | Verify API if version differs
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Data.CodonTable import TranslationError

# The silent trap: same sequence, two genetic codes, two different proteins.
# AGA is a stop in vertebrate mito (table 2) but Arg under the Standard code,
# and TGA is Trp in mito but a stop under Standard. Neither raises an error.
mito_gene = Seq('ATGAGATGA')
print('=== Same sequence, different genetic code ===')
print(f'Sequence:       {mito_gene}')
print(f'Standard (1):   {mito_gene.translate(table=1)}')
print(f'Vert Mito (2):  {mito_gene.translate(table=2)}')

# Table selection by integer id or registered name are equivalent
print('\n=== Table by id vs name ===')
print(mito_gene.translate(table=2))
print(mito_gene.translate(table='Vertebrate Mitochondrial'))

# Inspect a table: starts, stops, and a reassigned codon
print('\n=== Table 2 internals ===')
table = CodonTable.unambiguous_dna_by_id[2]
print(f'Start codons: {table.start_codons}')
print(f'Stop codons:  {table.stop_codons}')
print(f"TGA codes for: {table.forward_table['TGA']}")  # Trp, not stop

# cds=True validates a complete CDS and strips the terminal stop
print('\n=== cds=True on a valid CDS ===')
valid_cds = Seq('ATGTTTGGTTAA')
print(f'{valid_cds} -> {valid_cds.translate(cds=True)}')

# An alternative start (GTG) is read as M under the bacterial code
print('\n=== Alternative start under table 11 ===')
alt_start = Seq('GTGTTTGGTTAA')
print(f'{alt_start} -> {alt_start.translate(table=11, cds=True)}')

# cds=True turns silent defects into loud errors. TranslationError is the only
# expected failure here, so it is caught to demonstrate the message.
print('\n=== cds=True surfaces defects loudly ===')
for label, bad in [('no start', Seq('TTTGGTTAA')), ('length not x3', Seq('ATGTTTGG')), ('no stop', Seq('ATGTTTGGT'))]:
    try:
        bad.translate(cds=True)
    except TranslationError as e:
        print(f'{label}: {e}')
