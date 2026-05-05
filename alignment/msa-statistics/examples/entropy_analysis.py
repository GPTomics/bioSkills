'''Calculate Shannon entropy and Kullback-Leibler information content per column.'''
# Reference: biopython 1.83+ | Verify API if version differs

from Bio import AlignIO
from collections import Counter
import math

# Empirical amino acid background frequencies derived from Robinson & Robinson 1991 PNAS.
# Values match the commonly distributed set used in conservation-scoring code (e.g.
# Capra & Singh 2007 supplementary code). Sub-1% deviations from any single source
# table are not material for downstream KL-divergence ranking.
ROBINSON_BACKGROUND = {
    'A': 0.0780, 'R': 0.0512, 'N': 0.0427, 'D': 0.0530, 'C': 0.0193,
    'Q': 0.0419, 'E': 0.0629, 'G': 0.0738, 'H': 0.0224, 'I': 0.0526,
    'L': 0.0922, 'K': 0.0596, 'M': 0.0224, 'F': 0.0399, 'P': 0.0508,
    'S': 0.0712, 'T': 0.0584, 'W': 0.0133, 'Y': 0.0327, 'V': 0.0653,
}
DNA_UNIFORM = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}

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

def information_content(column, background, ignore_gaps=True):
    '''KL-divergence IC against a background distribution.

    Use ROBINSON_BACKGROUND for protein (Robinson & Robinson 1991 PNAS) or
    DNA_UNIFORM for nucleotide alignments.
    '''
    if ignore_gaps:
        column = column.replace('-', '')
    if not column:
        return 0.0
    counts = Counter(column)
    total = len(column)
    ic = 0.0
    for residue, count in counts.items():
        observed = count / total
        expected = background.get(residue, 1e-9)
        if observed > 0:
            ic += observed * math.log2(observed / expected)
    return ic

if __name__ == '__main__':
    alignment = AlignIO.read('alignment.fasta', 'fasta')
    print(f'Alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} columns\n')

    is_protein_guess = any(c in str(alignment[0].seq).upper() for c in 'EFILPQYW')
    background = ROBINSON_BACKGROUND if is_protein_guess else DNA_UNIFORM
    alphabet_label = 'protein (Robinson 1991 background)' if is_protein_guess else 'DNA (uniform background)'
    print(f'Treating as {alphabet_label}\n')

    print('Position   Entropy (bits)   IC (bits)')
    print('-' * 45)
    entropies = []
    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]
        ent = shannon_entropy(column)
        ic = information_content(column, background)
        entropies.append(ent)
        print(f'{i:5d}      {ent:8.3f}         {ic:8.3f}')

    avg_entropy = sum(entropies) / len(entropies)
    print(f'\nAverage entropy: {avg_entropy:.3f} bits')
    max_h = math.log2(20) if is_protein_guess else math.log2(4)
    print(f'Maximum possible entropy: {max_h:.3f} bits')
