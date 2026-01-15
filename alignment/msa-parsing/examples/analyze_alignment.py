'''Analyze alignment composition and conservation'''

from Bio import AlignIO
from collections import Counter

alignment = AlignIO.read('alignment.fasta', 'fasta')
print(f'Alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} columns\n')

print('Column composition (first 10 columns):')
for col_idx in range(min(10, alignment.get_alignment_length())):
    column = alignment[:, col_idx]
    counts = Counter(column)
    most_common = counts.most_common(1)[0]
    conservation = most_common[1] / len(alignment) * 100
    print(f'  Col {col_idx}: {dict(counts)} - {conservation:.0f}% conserved')
