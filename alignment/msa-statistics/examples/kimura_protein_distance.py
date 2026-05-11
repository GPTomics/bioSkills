'''Kimura empirical distance correction for protein alignments.

Kimura's protein distance corrects observed pairwise identity for multiple substitutions
at the same site: d = -log(1 - p - 0.2 * p^2). Valid up to ~85% divergence; saturates
beyond. Use Bio.Phylo.TreeConstruction.DistanceCalculator('blosum62') for matrix-based
alternatives or IQ-TREE for full ML distance estimation.
'''
# Reference: biopython 1.83+, numpy 1.26+ | Verify API if version differs

from Bio import AlignIO
import numpy as np
import math

def kimura_protein_distance(seq1, seq2):
    matches = sum(a == b and a != '-' and b != '-' for a, b in zip(seq1, seq2))
    aligned = sum(a != '-' and b != '-' for a, b in zip(seq1, seq2))
    if aligned == 0:
        return float('nan')
    p = 1 - matches / aligned
    if p >= 0.85:
        return float('inf')
    return -math.log(1 - p - 0.2 * p ** 2)

if __name__ == '__main__':
    alignment = AlignIO.read('alignment.fasta', 'fasta')
    n = len(alignment)
    distance = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = kimura_protein_distance(str(alignment[i].seq), str(alignment[j].seq))
            distance[i, j] = distance[j, i] = d

    print('Kimura protein distance matrix:')
    seq_ids = [r.id for r in alignment]
    for i, row in enumerate(distance):
        print(f'{seq_ids[i][:10]:>10}', ' '.join(f'{v:6.3f}' for v in row))

    mask = np.triu(np.ones_like(distance, dtype=bool), k=1)
    finite = distance[mask & np.isfinite(distance)]
    print(f'\nMean pairwise distance: {finite.mean():.3f}')
    print(f'Saturated pairs (distance >= 0.85 identity loss): {(~np.isfinite(distance) & mask).sum()}')
