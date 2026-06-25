#!/usr/bin/env python3
'''
RNA secondary structure prediction using the ViennaRNA Python API.
Treats the Boltzmann ensemble as the object: MFE, partition function,
centroid, MEA, per-base positional entropy, ensemble defect, and stochastic
sampling of alternative conformations, plus constrained and SHAPE-directed folding.
'''
# Reference: ViennaRNA 2.6+, matplotlib 3.8+, numpy 1.26+ | Verify API if version differs

import RNA
import numpy as np
import matplotlib.pyplot as plt


def fold_ensemble(sequence):
    '''Fold one sequence and return ensemble quantities, not just the MFE.'''
    fc = RNA.fold_compound(sequence)

    mfe_struct, mfe = fc.mfe()
    _, ensemble_g = fc.pf()                 # partition function; required before the quantities below
    centroid_struct, _ = fc.centroid()      # conservative, ensemble-representative
    mea_struct, _ = fc.MEA()                # maximum expected accuracy

    diversity = fc.mean_bp_distance()       # ensemble diversity: low = well-defined (read RELATIVE to length)
    defect = fc.ensemble_defect(mfe_struct) # expected wrongly-paired positions of the MFE vs the ensemble

    return {'mfe_structure': mfe_struct, 'mfe_energy': mfe, 'ensemble_energy': ensemble_g,
            'centroid_structure': centroid_struct, 'mea_structure': mea_struct,
            'ensemble_diversity': diversity, 'ensemble_defect': defect}


def positional_confidence(sequence):
    '''Per-base Shannon entropy over pairing partners; low entropy = confident position.'''
    fc = RNA.fold_compound(sequence)
    fc.pf()
    entropy = fc.positional_entropy()       # index 0 is a placeholder; positions are 1-indexed
    return list(entropy)[1:len(sequence) + 1]


def get_bpp_matrix(sequence):
    '''Base-pair probability matrix (the dot plot); pairs with bpp > 0.9 are trustworthy.'''
    fc = RNA.fold_compound(sequence)
    fc.pf()
    bpp = fc.bpp()
    n = len(sequence)
    matrix = np.zeros((n, n))
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            matrix[i - 1][j - 1] = bpp[i][j]
            matrix[j - 1][i - 1] = bpp[i][j]
    return matrix


def plot_bpp_dotplot(sequence, output_file='bpp_dotplot.png'):
    '''Render the base-pair probability matrix.'''
    matrix = get_bpp_matrix(sequence)
    n = len(sequence)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(matrix, cmap='YlOrRd', origin='lower', vmin=0, vmax=1)
    if n <= 80:
        ax.set_xticks(range(0, n, 10))
        ax.set_yticks(range(0, n, 10))
    ax.set_xlabel('Position')
    ax.set_ylabel('Position')
    ax.set_title('Base-pair probability matrix')
    plt.colorbar(ax.images[0], ax=ax, label='Probability')
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f'Saved dot plot to {output_file}')


def fold_with_constraints(sequence, forced_unpaired=None, forced_pairs=None):
    '''Fold with hard positional constraints (1-indexed).'''
    fc = RNA.fold_compound(sequence)
    for pos in forced_unpaired or []:
        fc.hc_add_up(pos, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
    for i, j in forced_pairs or []:
        fc.hc_add_bp(i, j, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
    return fc.mfe()


def fold_with_shape(sequence, reactivities, m=1.8, b=-0.6):
    '''
    Fold under a Deigan SHAPE pseudo-energy restraint.

    reactivities is 1-indexed: index 0 is a -999 placeholder; -999 elsewhere means no data, not zero.
    m=1.8, b=-0.6 is the standard SHAPE pair (Hajdin 2013, the ViennaRNA default).
    '''
    fc = RNA.fold_compound(sequence)
    fc.sc_add_SHAPE_deigan(reactivities, m, b)
    return fc.mfe()


if __name__ == '__main__':
    trna = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'

    print('=== Ensemble folding ===')
    for key, val in fold_ensemble(trna).items():
        print(f'  {key}: {round(val, 2) if isinstance(val, float) else val}')

    print('\n=== Per-base confidence (positional entropy) ===')
    entropy = positional_confidence(trna)
    print(f'  mean entropy: {np.mean(entropy):.3f}  (low = confident; high = ambiguous)')

    print('\n=== Constrained folding (force the acceptor stem open) ===')
    unconstrained = fold_ensemble(trna)['mfe_structure']
    constrained, c_mfe = fold_with_constraints(trna, forced_unpaired=[2, 3, 4])
    print(f'  unconstrained: {unconstrained}')
    print(f'  constrained:   {constrained} ({c_mfe:.2f} kcal/mol)')
    print(f'  base-pair distance: {RNA.bp_distance(unconstrained, constrained)}')

    print('\n=== SHAPE-directed folding ===')
    # Toy reactivities: low where the unconstrained fold pairs, high where it does not.
    react = [-999.0] + [0.15 if c in '()' else 0.8 for c in unconstrained]
    shape_struct, shape_mfe = fold_with_shape(trna, react)
    # The reported energy INCLUDES the SHAPE pseudo-energy, so it is NOT comparable to the unconstrained MFE.
    print(f'  SHAPE-guided: {shape_struct} ({shape_mfe:.2f} kcal/mol, includes SHAPE pseudo-energy)')

    print('\n=== Base-pair probability dot plot ===')
    plot_bpp_dotplot(trna, 'trna_bpp_dotplot.png')

    print('\nDone.')
