'''Compute Henikoff position-based sequence weights for an MSA.

Reference: Henikoff S, Henikoff JG. 1994. Position-based sequence weights. JMB 243:574-578.
Each column contributes 1 / (k * n) per sequence, where k is the number of distinct residues
and n is the count of the sequence's residue at that column.

Gap handling: columns containing any gap character are skipped entirely. This matches HMMER's
esl-weight convention of restricting position-based weighting to ungapped (consensus) columns.
The alternative (treating '-' as a 21st residue) inflates weights toward sequences with unique
gap patterns. For MSAs where many sequences are fragmentary, prefer pyhmmer's compute_weights
or restrict to a match-state-only sub-alignment first.
'''
# Reference: biopython 1.83+, numpy 1.26+ | Verify API if version differs

from Bio import AlignIO
import numpy as np

def henikoff_weights(alignment):
    n_seqs = len(alignment)
    n_cols = alignment.get_alignment_length()
    seq_array = np.array([list(str(r.seq)) for r in alignment])
    weights = np.zeros(n_seqs)
    for col_idx in range(n_cols):
        column = seq_array[:, col_idx]
        residues, inverse, counts = np.unique(column, return_inverse=True, return_counts=True)
        if '-' in residues:
            continue
        n_distinct = len(residues)
        per_seq_weight = 1.0 / (n_distinct * counts[inverse])
        weights += per_seq_weight
    return weights / weights.sum()

if __name__ == '__main__':
    alignment = AlignIO.read('alignment.fasta', 'fasta')
    weights = henikoff_weights(alignment)
    for record, weight in zip(alignment, weights):
        print(f'{record.id}: {weight:.4f}')
    print(f'\nTotal weight: {weights.sum():.4f}  Effective sequences: {1 / (weights ** 2).sum():.2f}')
