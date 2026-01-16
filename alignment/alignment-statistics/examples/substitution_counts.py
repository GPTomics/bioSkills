'''Count substitutions observed in alignment'''

from Bio import AlignIO
from collections import defaultdict

def substitution_counts(alignment):
    counts = defaultdict(int)
    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx]
        chars = [c for c in column if c != '-']
        for i, c1 in enumerate(chars):
            for c2 in chars[i+1:]:
                if c1 != c2:
                    pair = tuple(sorted([c1, c2]))
                    counts[pair] += 1
    return dict(counts)

alignment = AlignIO.read('alignment.fasta', 'fasta')
print(f'Alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} columns\n')

subs = substitution_counts(alignment)
total_subs = sum(subs.values())

print(f'Total substitutions observed: {total_subs}\n')
print('Substitution counts (sorted by frequency):')
for pair, count in sorted(subs.items(), key=lambda x: -x[1]):
    pct = count / total_subs * 100
    print(f'  {pair[0]} <-> {pair[1]}: {count:5d} ({pct:5.1f}%)')

transitions = sum(v for k, v in subs.items() if set(k) in [{'A', 'G'}, {'C', 'T'}])
transversions = total_subs - transitions
print(f'\nTransitions: {transitions}')
print(f'Transversions: {transversions}')
if transversions > 0:
    print(f'Ti/Tv ratio: {transitions/transversions:.2f}')
