'''Basic examples of creating and using Seq objects, plus post-1.78 traps'''
# Reference: biopython 1.83+ | Verify API if version differs
from Bio.Seq import Seq, MutableSeq, UndefinedSequenceError
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

# Trap: Seq is immutable since 1.78 (use MutableSeq to edit in place)
print('\n=== Immutability Trap ===')
try:
    dna[0] = 'C'
except TypeError as err:
    print(f'Seq edit blocked: {err}')

# Trap: Seq is bytes-backed since 1.79, NOT a str subclass
print('\n=== isinstance Trap ===')
print(f'isinstance(dna, str): {isinstance(dna, str)}')
print(f'isinstance(dna, (Seq, MutableSeq)): {isinstance(dna, (Seq, MutableSeq))}')

# Trap: undefined sequence has a length but no readable content (lazy parsers)
print('\n=== Undefined Sequence Trap ===')
undef = Seq(None, length=20)
print(f'Length known: {len(undef)}')
try:
    str(undef)
except UndefinedSequenceError as err:
    print(f'Content read blocked: {err}')
