'''Choosing an access pattern: streaming, in-memory, and indexed reads'''
# Reference: biopython 1.83+ | Verify API if version differs
import os
from Bio import SeqIO

# Streaming: parse() is a one-pass generator that exhausts after the first loop
print('=== parse() is one-pass ===')
records = SeqIO.parse('sample.fasta', 'fasta')
first_pass = [r.id for r in records]
second_pass = [r.id for r in records]  # generator already exhausted, no error
print(f'First pass:  {first_pass}')
print(f'Second pass: {second_pass}  (silently empty)')

# In-memory random access for small files
print('\n=== to_dict() for small files ===')
by_id = SeqIO.to_dict(SeqIO.parse('sample.fasta', 'fasta'))
print(f'seq2 length: {len(by_id["seq2"].seq)} bp')

# Indexed random access without loading everything (re-parses on each access)
print('\n=== index() for large files ===')
index = SeqIO.index('sample.fasta', 'fasta')
print(f'seq3 length: {len(index["seq3"].seq)} bp')
index.close()

# Persistent on-disk index that survives across sessions
print('\n=== index_db() persists across runs ===')
db_path = 'demo_index.sqlite'
db = SeqIO.index_db(db_path, 'sample.fasta', 'fasta')
print(f'Indexed IDs: {list(db.keys())}')
db.close()
os.remove(db_path)

# id is the first header token; description is the whole header
print('\n=== id vs name vs description ===')
record = next(SeqIO.parse('sample.fasta', 'fasta'))
print(f'id:          {record.id}')
print(f'name:        {record.name}')
print(f'description: {record.description}')
