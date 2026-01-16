'''Splitting sequences into chunks and windows'''
from Bio.Seq import Seq

seq = Seq('ATGCGATCGATCGATCGATCGATCG')
print(f'Original: {seq} ({len(seq)} bp)')

# Split into codons
print('\n=== Split into Codons ===')
def split_codons(seq):
    return [seq[i:i+3] for i in range(0, len(seq) - len(seq) % 3, 3)]

codons = split_codons(seq)
for i, codon in enumerate(codons):
    print(f'Codon {i+1}: {codon}')

# Split into fixed chunks
print('\n=== Fixed-Size Chunks (10 bp) ===')
def chunk_sequence(seq, size):
    return [seq[i:i+size] for i in range(0, len(seq), size)]

chunks = chunk_sequence(seq, 10)
for i, chunk in enumerate(chunks):
    print(f'Chunk {i+1}: {chunk} ({len(chunk)} bp)')

# Sliding windows
print('\n=== Sliding Windows (6 bp, step 3) ===')
def sliding_windows(seq, window_size, step=1):
    for i in range(0, len(seq) - window_size + 1, step):
        yield i, seq[i:i + window_size]

for pos, window in sliding_windows(seq, 6, 3):
    print(f'Position {pos}: {window}')

# Get flanking regions
print('\n=== Flanking Regions ===')
def get_flanking(seq, position, flank_size):
    start = max(0, position - flank_size)
    end = min(len(seq), position + flank_size + 1)
    return start, end, seq[start:end]

pos = 12
flank = 5
start, end, region = get_flanking(seq, pos, flank)
print(f'Position {pos} with {flank} bp flanks: positions {start}-{end}')
print(f'Sequence: {region}')
