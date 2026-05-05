'''Compute effective sequence number (Neff) by similarity-based clustering.

Cluster sequences sharing >= identity_threshold of non-gap positions and assign each
sequence the inverse of its cluster size as a weight; Neff is the sum of weights.
Used by HMMER and AlphaFold as an MSA-depth metric.
'''
# Reference: biopython 1.83+, numpy 1.26+ | Verify API if version differs

from Bio import AlignIO
import numpy as np

def neff(alignment, identity_threshold=0.62):
    seq_array = np.array([list(str(r.seq)) for r in alignment])
    n = len(alignment)
    cluster_size = np.ones(n)
    for i in range(n):
        for j in range(i + 1, n):
            mask = (seq_array[i] != '-') & (seq_array[j] != '-')
            length = mask.sum()
            if not length:
                continue
            shared = ((seq_array[i] == seq_array[j]) & mask).sum()
            if shared / length >= identity_threshold:
                cluster_size[i] += 1
                cluster_size[j] += 1
    weights = 1.0 / cluster_size
    return weights.sum()

if __name__ == '__main__':
    alignment = AlignIO.read('alignment.fasta', 'fasta')
    n_eff_protein = neff(alignment, identity_threshold=0.62)
    n_eff_nucleotide = neff(alignment, identity_threshold=0.80)
    length = alignment.get_alignment_length()

    print(f'Sequences: {len(alignment)}')
    print(f'Length: {length}')
    print(f'Neff (62% threshold, protein convention): {n_eff_protein:.2f}')
    print(f'Neff/L: {n_eff_protein / length:.3f}')
    print(f'Neff (80% threshold, nucleotide convention): {n_eff_nucleotide:.2f}')
    print('Rule of thumb: Neff/L > 0.5 sufficient for direct-coupling-analysis contact prediction.')
