'''GC content analysis across sequences and windows'''
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from collections import Counter

# GC content of a single sequence
print('=== Single Sequence ===')
seq = Seq('ATGCGATCGATCGATCGATCGATGCGCGCGCGCATATATATATAT')
print(f'Sequence: {seq}')
print(f'GC content: {gc_fraction(seq) * 100:.1f}%')

# Sliding window GC
print('\n=== Sliding Window GC (10 bp, step 5) ===')
def gc_windows(seq, window_size=10, step=5):
    results = []
    for i in range(0, len(seq) - window_size + 1, step):
        window = seq[i:i + window_size]
        gc = gc_fraction(window) * 100
        results.append((i, i + window_size, gc))
    return results

for start, end, gc in gc_windows(seq, 10, 5):
    bar = '#' * int(gc / 5)
    print(f'{start:3d}-{end:3d}: {gc:5.1f}% {bar}')

# Codon usage
print('\n=== Codon Usage ===')
coding_seq = Seq('ATGCGATCGATCGATCGATCGTAA')

def codon_usage(seq):
    seq_str = str(seq)
    codons = [seq_str[i:i+3] for i in range(0, len(seq_str) - 2, 3)]
    counts = Counter(codons)
    total = sum(counts.values())
    return {codon: (count, count / total * 100) for codon, count in counts.items()}

usage = codon_usage(coding_seq)
print(f'Coding sequence: {coding_seq}')
for codon, (count, pct) in sorted(usage.items()):
    print(f'{codon}: {count} ({pct:.1f}%)')

# Dinucleotide frequencies
print('\n=== Dinucleotide Frequencies ===')
def dinucleotide_freq(seq):
    seq_str = str(seq)
    dinucs = [seq_str[i:i+2] for i in range(len(seq_str) - 1)]
    counts = Counter(dinucs)
    total = sum(counts.values())
    return {di: count / total * 100 for di, count in counts.items()}

freq = dinucleotide_freq(seq)
for di, pct in sorted(freq.items(), key=lambda x: x[1], reverse=True)[:8]:
    print(f'{di}: {pct:.1f}%')

# CpG ratio
print('\n=== CpG Analysis ===')
cpg_count = str(seq).count('CG')
expected_cpg = (str(seq).count('C') * str(seq).count('G')) / len(seq)
cpg_ratio = cpg_count / expected_cpg if expected_cpg > 0 else 0
print(f'CpG observed: {cpg_count}')
print(f'CpG expected: {expected_cpg:.1f}')
print(f'CpG O/E ratio: {cpg_ratio:.2f}')
