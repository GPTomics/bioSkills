'''Fetch sequences from NCBI with Bio.Entrez and Bio.SeqIO'''
from Bio import Entrez, SeqIO

Entrez.email = 'your.email@example.com'

# Fetch single FASTA
print('=== Fetch FASTA ===')
handle = Entrez.efetch(db='nucleotide', id='NM_007294', rettype='fasta', retmode='text')
record = SeqIO.read(handle, 'fasta')
handle.close()
print(f"{record.id}: {len(record.seq)} bp")
print(f"Sequence: {record.seq[:60]}...")

# Fetch GenBank with features
print('\n=== Fetch GenBank ===')
handle = Entrez.efetch(db='nucleotide', id='NM_007294', rettype='gb', retmode='text')
record = SeqIO.read(handle, 'genbank')
handle.close()
print(f"Description: {record.description}")
print(f"Features: {len(record.features)}")
for f in record.features[:5]:
    print(f"  {f.type}: {f.location}")

# Fetch multiple sequences
print('\n=== Fetch Multiple ===')
ids = ['NM_007294', 'NM_000059', 'NM_000546']
handle = Entrez.efetch(db='nucleotide', id=','.join(ids), rettype='fasta', retmode='text')
records = list(SeqIO.parse(handle, 'fasta'))
handle.close()
for r in records:
    print(f"{r.id}: {len(r.seq)} bp")

# Save to file
print('\n=== Save to File ===')
handle = Entrez.efetch(db='nucleotide', id=','.join(ids), rettype='fasta', retmode='text')
with open('downloaded.fasta', 'w') as out:
    out.write(handle.read())
handle.close()
print('Saved to downloaded.fasta')
