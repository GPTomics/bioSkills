'''Calculate Shannon entropy and information content'''

from Bio import AlignIO
from collections import Counter
import math

def shannon_entropy(column, ignore_gaps=True):
    if ignore_gaps:
        column = column.replace('-', '')
    if not column:
        return 0.0
    counts = Counter(column)
    total = len(column)
    entropy = 0.0
    for count in counts.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy

def information_content(column, alphabet_size=4, ignore_gaps=True):
    max_entropy = math.log2(alphabet_size)
    observed_entropy = shannon_entropy(column, ignore_gaps)
    return max_entropy - observed_entropy

alignment = AlignIO.read('alignment.fasta', 'fasta')
print(f'Alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} columns\n')

print('Position   Entropy (bits)   Info Content')
print('-' * 45)
for i in range(alignment.get_alignment_length()):
    column = alignment[:, i]
    ent = shannon_entropy(column)
    ic = information_content(column, alphabet_size=4)
    print(f'{i:5d}      {ent:8.3f}         {ic:8.3f}')

entropies = [shannon_entropy(alignment[:, i]) for i in range(alignment.get_alignment_length())]
avg_entropy = sum(entropies) / len(entropies)
print(f'\nAverage entropy: {avg_entropy:.3f} bits')
print(f'Max possible entropy (DNA): {math.log2(4):.3f} bits')
