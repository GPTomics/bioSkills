'''Basic examples of creating and using Seq objects'''
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord

# Create Seq from string
print('=== Creating Seq Objects ===')
dna = Seq('ATGCGATCGATCG')
print(f'DNA: {dna}')
print(f'Length: {len(dna)}')
print(f'First 3 bases: {dna[:3]}')

# String-like operations
print('\n=== String Operations ===')
print(f'Count G: {dna.count("G")}')
print(f'Find ATG: position {dna.find("ATG")}')
print(f'Contains GAT: {"GAT" in dna}')

# MutableSeq for in-place modifications
print('\n=== MutableSeq ===')
mut_seq = MutableSeq('ATGCGATCG')
print(f'Original: {mut_seq}')
mut_seq[0] = 'G'
print(f'After mut_seq[0] = "G": {mut_seq}')
mut_seq.extend('TTT')
print(f'After extend: {mut_seq}')

# Convert back to immutable
final_seq = Seq(mut_seq)
print(f'Final Seq: {final_seq}')

# Create SeqRecord
print('\n=== SeqRecord ===')
record = SeqRecord(Seq('ATGCGATCGATCG'), id='gene1', description='Example gene sequence')
print(f'ID: {record.id}')
print(f'Description: {record.description}')
print(f'Sequence: {record.seq}')
