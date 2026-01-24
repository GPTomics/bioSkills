'''Calculate pairwise identity matrix for alignment'''

from Bio import AlignIO
import numpy as np

def pairwise_identity(seq1, seq2):
    matches = sum(a == b and a != '-' for a, b in zip(seq1, seq2))
    aligned_positions = sum(a != '-' or b != '-' for a, b in zip(seq1, seq2))
    return matches / aligned_positions if aligned_positions > 0 else 0

alignment = AlignIO.read('alignment.fasta', 'fasta')
n = len(alignment)
seq_ids = [r.id for r in alignment]

matrix = np.zeros((n, n))
for i in range(n):
    for j in range(i, n):
        seq_i = str(alignment[i].seq)
        seq_j = str(alignment[j].seq)
        ident = pairwise_identity(seq_i, seq_j)
        matrix[i, j] = matrix[j, i] = ident

print('Pairwise Identity Matrix (%):')
print(f'{"":>12}', ' '.join(f'{s[:8]:>8}' for s in seq_ids))
for i, row in enumerate(matrix):
    print(f'{seq_ids[i][:12]:>12}', ' '.join(f'{v*100:>7.1f}%' for v in row))

avg_identity = (matrix.sum() - n) / (n * (n - 1))
print(f'\nAverage pairwise identity: {avg_identity*100:.1f}%')
