'''Basic sequence writing examples'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Create simple records
records = [
    SeqRecord(Seq('ATGCGATCGATCGATCGATCG'), id='seq1', description='First sequence'),
    SeqRecord(Seq('GCTAGCTAGCTAGCTAGCTA'), id='seq2', description='Second sequence'),
    SeqRecord(Seq('TTAATTAATTAATTAATTAA'), id='seq3', description='Third sequence')
]

# Write to FASTA
count = SeqIO.write(records, 'output.fasta', 'fasta')
print(f'Wrote {count} records to output.fasta')

# Get formatted string without writing
for record in records:
    print(record.format('fasta'))
