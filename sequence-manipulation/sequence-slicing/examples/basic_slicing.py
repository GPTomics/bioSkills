'''Basic sequence slicing and indexing operations'''
from Bio.Seq import Seq

seq = Seq('ATGCGATCGATCGATCGATCG')
print(f'Original sequence: {seq}')
print(f'Length: {len(seq)}')

# Indexing single positions
print('\n=== Indexing ===')
print(f'First base (seq[0]): {seq[0]}')
print(f'Last base (seq[-1]): {seq[-1]}')
print(f'Fifth base (seq[4]): {seq[4]}')

# Basic slicing
print('\n=== Slicing ===')
print(f'First 6 bases (seq[:6]): {seq[:6]}')
print(f'Last 6 bases (seq[-6:]): {seq[-6:]}')
print(f'Positions 3-9 (seq[3:9]): {seq[3:9]}')

# Step slicing
print('\n=== Step Slicing ===')
print(f'Every 3rd base (seq[::3]): {seq[::3]}')
print(f'Reversed (seq[::-1]): {seq[::-1]}')

# Slicing returns Seq objects
print('\n=== Return Types ===')
sliced = seq[0:6]
print(f'Type of slice: {type(sliced)}')
print(f'Can translate: {sliced.translate()}')

# 1-based coordinate conversion
print('\n=== 1-Based Coordinates ===')
def extract_1based(seq, start, end):
    '''Extract using 1-based inclusive coordinates'''
    return seq[start - 1:end]

region = extract_1based(seq, 1, 6)  # Positions 1-6 in biology notation
print(f'Positions 1-6 (1-based): {region}')
