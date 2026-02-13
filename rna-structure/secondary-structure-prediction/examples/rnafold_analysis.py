#!/usr/bin/env python3
'''
RNA secondary structure prediction using ViennaRNA Python API.
Demonstrates MFE folding, partition function analysis, constrained folding,
and structure comparison.
'''
# Reference: infernal 1.1+, matplotlib 3.8+, numpy 1.26+ | Verify API if version differs

import RNA
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def fold_sequence(sequence):
    '''Fold a single RNA sequence and return all structure predictions.'''
    fc = RNA.fold_compound(sequence)

    mfe_struct, mfe = fc.mfe()
    pf_struct, pf_energy = fc.pf()
    centroid_struct, centroid_dist = fc.centroid()
    mea_struct, mea_val = fc.MEA()

    # Ensemble diversity: low values indicate well-defined structure
    ensemble_diversity = fc.mean_bp_distance()

    return {
        'mfe_structure': mfe_struct,
        'mfe_energy': mfe,
        'pf_structure': pf_struct,
        'ensemble_energy': pf_energy,
        'centroid_structure': centroid_struct,
        'centroid_distance': centroid_dist,
        'mea_structure': mea_struct,
        'mea_value': mea_val,
        'ensemble_diversity': ensemble_diversity
    }


def get_bpp_matrix(sequence):
    '''Compute base-pair probability matrix.'''
    fc = RNA.fold_compound(sequence)
    fc.pf()
    bpp = fc.bpp()

    n = len(sequence)
    matrix = np.zeros((n, n))
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            matrix[i-1][j-1] = bpp[i][j]
            matrix[j-1][i-1] = bpp[i][j]
    return matrix


def plot_bpp_dotplot(sequence, output_file='bpp_dotplot.png'):
    '''Generate base-pair probability dot plot.'''
    matrix = get_bpp_matrix(sequence)
    n = len(sequence)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(matrix, cmap='YlOrRd', origin='lower', vmin=0, vmax=1)

    if n <= 80:
        tick_positions = range(0, n, 10)
        ax.set_xticks(tick_positions)
        ax.set_yticks(tick_positions)

    ax.set_xlabel('Position')
    ax.set_ylabel('Position')
    ax.set_title('Base-pair probability matrix')
    plt.colorbar(ax.images[0], ax=ax, label='Probability')
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f'Saved dot plot to {output_file}')


def fold_with_constraints(sequence, forced_unpaired=None, forced_pairs=None):
    '''
    Fold with positional constraints.

    forced_unpaired: list of 1-indexed positions to force unpaired
    forced_pairs: list of (i, j) 1-indexed position tuples to force paired
    '''
    md = RNA.md()
    md.uniq_ML = 1
    fc = RNA.fold_compound(sequence, md)

    if forced_unpaired:
        for pos in forced_unpaired:
            fc.hc_add_up(pos, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)

    if forced_pairs:
        for i, j in forced_pairs:
            fc.hc_add_bp(i, j, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)

    structure, mfe = fc.mfe()
    return structure, mfe


def fold_with_shape(sequence, shape_reactivities, m=1.8, b=-0.6):
    '''
    Fold with SHAPE reactivity constraints.

    shape_reactivities: list of per-nucleotide reactivities (1-indexed).
                        Use -999 for missing data positions.
    m, b: Deigan et al. (2009) slope and intercept parameters.
          m=1.8, b=-0.6 are standard for SHAPE-MaP.
    '''
    fc = RNA.fold_compound(sequence)
    fc.sc_add_SHAPE_deigan(shape_reactivities, m, b)
    structure, mfe = fc.mfe()
    return structure, mfe


def compare_structures(struct1, struct2):
    '''Compare two structures using base-pair and tree edit distances.'''
    bp_dist = RNA.bp_distance(struct1, struct2)

    tree1 = RNA.make_tree(RNA.expand_Full(struct1))
    tree2 = RNA.make_tree(RNA.expand_Full(struct2))
    tree_dist = RNA.tree_edit_distance(tree1, tree2)

    return {'bp_distance': bp_dist, 'tree_edit_distance': tree_dist}


def structure_to_pairs(structure):
    '''Convert dot-bracket to list of base pairs (1-indexed).'''
    pt = RNA.ptable(structure)
    return [(i, pt[i]) for i in range(1, len(pt)) if pt[i] > i]


def compute_mfe_zscore(sequence, n_shuffles=100):
    '''
    Compute z-score of MFE relative to dinucleotide-shuffled controls.

    z-score < -2.0: sequence folds significantly better than random,
    suggesting functional structure.
    '''
    _, native_mfe = RNA.fold(sequence)

    shuffled_mfes = []
    for _ in range(n_shuffles):
        shuffled = RNA.sequence_shuffle(sequence)
        _, shuf_mfe = RNA.fold(shuffled)
        shuffled_mfes.append(shuf_mfe)

    mean_shuf = np.mean(shuffled_mfes)
    std_shuf = np.std(shuffled_mfes)

    # Avoid division by zero for very short sequences
    zscore = (native_mfe - mean_shuf) / std_shuf if std_shuf > 0 else 0.0
    return zscore, native_mfe, mean_shuf, std_shuf


if __name__ == '__main__':
    # Example: tRNA-Phe sequence
    trna_seq = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'

    print('=== MFE Folding ===')
    results = fold_sequence(trna_seq)
    for key, val in results.items():
        if isinstance(val, float):
            print(f'  {key}: {val:.2f}')
        else:
            print(f'  {key}: {val}')

    print('\n=== Structure Comparison ===')
    struct_unconstrained = results['mfe_structure']
    struct_constrained, constrained_mfe = fold_with_constraints(
        trna_seq, forced_unpaired=[35, 36, 37]
    )
    print(f'  Unconstrained: {struct_unconstrained}')
    print(f'  Constrained:   {struct_constrained} ({constrained_mfe:.2f} kcal/mol)')

    distances = compare_structures(struct_unconstrained, struct_constrained)
    print(f'  BP distance: {distances["bp_distance"]}')
    print(f'  Tree edit distance: {distances["tree_edit_distance"]}')

    print('\n=== MFE Z-score ===')
    zscore, native, mean_shuf, std_shuf = compute_mfe_zscore(trna_seq, n_shuffles=50)
    print(f'  Native MFE: {native:.2f} kcal/mol')
    print(f'  Shuffled mean: {mean_shuf:.2f} +/- {std_shuf:.2f}')
    print(f'  Z-score: {zscore:.2f}')

    print('\n=== Base-pair Probability Dot Plot ===')
    plot_bpp_dotplot(trna_seq, 'trna_bpp_dotplot.png')

    print('\nDone.')
