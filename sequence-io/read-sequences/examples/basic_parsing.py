'''Basic examples of reading sequences with Bio.SeqIO'''
from Bio import SeqIO

# Parse multiple sequences from FASTA
print('=== Parsing FASTA ===')
for record in SeqIO.parse('sample.fasta', 'fasta'):
    print(f'{record.id}: {len(record.seq)} bp')

# Read single sequence
print('\n=== Reading Single Sequence ===')
single = SeqIO.read('single.fasta', 'fasta')
print(f'{single.id}: {single.seq[:20]}...')

# Load all into list
print('\n=== Loading All Sequences ===')
all_records = list(SeqIO.parse('sample.fasta', 'fasta'))
print(f'Loaded {len(all_records)} sequences')

# Count without loading
print('\n=== Counting Sequences ===')
count = sum(1 for _ in SeqIO.parse('sample.fasta', 'fasta'))
print(f'File contains {count} sequences')

# Get just the IDs
print('\n=== Extracting IDs ===')
ids = [r.id for r in SeqIO.parse('sample.fasta', 'fasta')]
print(f'IDs: {ids}')
