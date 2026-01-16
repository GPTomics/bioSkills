'''Basic pattern searching in sequences'''
from Bio.Seq import Seq

seq = Seq('ATGCGAATTCGATCGAATTCGATCGGATCC')
print(f'Sequence: {seq}')
print(f'Length: {len(seq)} bp')

# Find first occurrence
print('\n=== find() - First Occurrence ===')
pos = seq.find('GAATTC')
print(f'First GAATTC at position: {pos}')

# Count occurrences
print('\n=== count() - Total Count ===')
n = seq.count('GAATTC')
print(f'GAATTC occurs {n} times')

# Find all occurrences
print('\n=== Find All Positions ===')
def find_all(seq, pattern):
    positions = []
    seq_str, pattern = str(seq), str(pattern)
    pos = seq_str.find(pattern)
    while pos != -1:
        positions.append(pos)
        pos = seq_str.find(pattern, pos + 1)
    return positions

positions = find_all(seq, 'GAATTC')
print(f'GAATTC positions: {positions}')

# Extract matches with context
print('\n=== Matches with Context ===')
pattern = 'GAATTC'
for pos in positions:
    start = max(0, pos - 3)
    end = min(len(seq), pos + len(pattern) + 3)
    context = seq[start:end]
    print(f'Position {pos}: ...{context}...')

# Search for multiple patterns
print('\n=== Multiple Patterns ===')
patterns = ['GAATTC', 'GGATCC', 'AAGCTT']
for pat in patterns:
    pos_list = find_all(seq, pat)
    status = f'at {pos_list}' if pos_list else 'not found'
    print(f'{pat}: {status}')
