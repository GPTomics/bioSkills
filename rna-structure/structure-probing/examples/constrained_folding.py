#!/usr/bin/env python3
'''
SHAPE/DMS reactivity-restrained RNA folding.
Loads a ShapeMapper2 profile, folds with and without the reactivity restraint,
and reports how the restraint changes the structure. Reactivity is a SOFT restraint
on folding, not a structure, and reports flexibility, not pairing per se.
'''
# Reference: ViennaRNA 2.6+, matplotlib 3.8+, numpy 1.26+, pandas 2.2+ | Verify API if version differs

import RNA
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def load_shape_profile(profile_file):
    '''
    Load a ShapeMapper2 <name>_<RNA>_profile.txt.

    The base character is the 'Sequence' column ('Nucleotide' is the 1-based position integer).
    Fold with 'Norm_profile' (normalized), not the raw 'Reactivity_profile'.
    '''
    df = pd.read_csv(profile_file, sep='\t')
    sequence = ''.join(df['Sequence'].tolist())
    react_col = 'Norm_profile' if 'Norm_profile' in df.columns else 'Reactivity_profile'
    reactivities = [r if pd.notna(r) else -999 for r in df[react_col].tolist()]
    return sequence, reactivities


def fold_unconstrained(sequence):
    '''Fold without experimental restraints.'''
    fc = RNA.fold_compound(sequence)
    structure, mfe = fc.mfe()
    fc.pf()
    centroid, _ = fc.centroid()
    return {'mfe': structure, 'mfe_energy': mfe, 'centroid': centroid, 'diversity': fc.mean_bp_distance()}


def fold_shape_constrained(sequence, reactivities, m=1.8, b=-0.6):
    '''
    Fold under a Deigan SHAPE pseudo-energy restraint.

    m=1.8, b=-0.6 is the standard SHAPE pair (Hajdin et al. 2013, the ViennaRNA default),
    NOT Deigan et al. 2009's own values (m=2.6, b=-0.8). There is no separate DMS-specific
    standard pair; reuse 1.8/-0.6 (Cordero et al. 2012 showed SHAPE parameters transfer) or tune,
    and mask G/U to -999 (DMS carries no Watson-Crick signal there). The vector is 1-indexed:
    index 0 is a -999 placeholder; -999 elsewhere means no data, not zero reactivity.
    '''
    shape_data = [-999.0] + [r if r != -999 else -999.0 for r in reactivities]
    fc = RNA.fold_compound(sequence)
    fc.sc_add_SHAPE_deigan(shape_data, m, b)
    structure, mfe = fc.mfe()
    fc.pf()
    centroid, _ = fc.centroid()
    return {'mfe': structure, 'mfe_energy': mfe, 'centroid': centroid, 'diversity': fc.mean_bp_distance()}


def reactivity_structure_agreement(structure, reactivities, low=0.4, high=0.85):
    '''
    Heuristic agreement: paired positions tend to have low reactivity, unpaired high.
    Defaults are the standard Weeks-lab SHAPE bins (unreactive <0.4, intermediate 0.4-0.85,
    reactive >0.85). This is a sanity check, not validation -- low reactivity can also reflect
    tertiary contacts or protein/ligand protection, not pairing.
    '''
    agree, total = 0, 0
    for char, react in zip(structure, reactivities):
        if react == -999 or low <= react <= high:
            continue
        total += 1
        paired = char in '()'
        if (paired and react < low) or (not paired and react > high):
            agree += 1
    return agree / total if total > 0 else 0.0


def compare_predictions(sequence, reactivities):
    '''Compare unrestrained and reactivity-restrained folds.'''
    unconstrained = fold_unconstrained(sequence)
    constrained = fold_shape_constrained(sequence, reactivities)
    bp_dist = RNA.bp_distance(unconstrained['mfe'], constrained['mfe'])

    print('=== Unrestrained ===')
    print(f'  MFE: {unconstrained["mfe"]} ({unconstrained["mfe_energy"]:.2f} kcal/mol)')
    print(f'  agreement with reactivities: {reactivity_structure_agreement(unconstrained["mfe"], reactivities):.1%}')
    print('\n=== Reactivity-restrained ===')
    # The restrained energy INCLUDES the SHAPE pseudo-energy, so it is NOT comparable to the unrestrained MFE.
    print(f'  MFE: {constrained["mfe"]} ({constrained["mfe_energy"]:.2f} kcal/mol, includes SHAPE pseudo-energy)')
    print(f'  agreement with reactivities: {reactivity_structure_agreement(constrained["mfe"], reactivities):.1%}')
    print(f'\nbase-pair distance between folds: {bp_dist} (compare structures, not the two energies)')
    return unconstrained, constrained


def plot_reactivity(sequence, structure, reactivities, output_file='reactivity_structure.png'):
    '''Bar plot of reactivity colored by predicted pairing status.'''
    n = len(sequence)
    positions = np.arange(1, n + 1)
    valid = np.array([r != -999 for r in reactivities])
    react_arr = np.array([r if r != -999 else 0 for r in reactivities])
    paired = np.array([c in '()' for c in structure])
    colors = np.where(paired, '#4169E1', '#FF4500')

    fig, ax = plt.subplots(figsize=(max(8, n * 0.08), 4))
    ax.bar(positions[valid], react_arr[valid], color=colors[valid], width=1.0, edgecolor='none')
    # Standard SHAPE reactivity bins (Weeks lab): unreactive <0.4, intermediate 0.4-0.85, reactive >0.85.
    # For the conventional figure (reactivity-colored arc diagram on the structure) use ShapeMapper arcPlot or forna.
    ax.axhline(0.4, color='gray', linestyle='--', linewidth=0.5)
    ax.axhline(0.85, color='gray', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Position')
    ax.set_ylabel('Normalized reactivity')
    ax.set_title('Reactivity (blue = predicted paired, red = predicted unpaired)')
    ax.set_xlim(0, n + 1)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f'Saved to {output_file}')


if __name__ == '__main__':
    profile_file = 'results/my_rna_my_rna_profile.txt'
    try:
        sequence, reactivities = load_shape_profile(profile_file)
    except FileNotFoundError:
        print(f'Profile not found: {profile_file}; using example tRNA reactivities.')
        sequence = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'
        rng = np.random.default_rng(42)
        reactivities = []
        for i in range(len(sequence)):
            if i < 7 or (25 < i < 43) or i > 65:
                reactivities.append(round(float(rng.uniform(0.0, 0.3)), 3))
            elif (14 <= i <= 21) or (32 <= i <= 37):
                reactivities.append(round(float(rng.uniform(0.5, 1.2)), 3))
            else:
                reactivities.append(round(float(rng.uniform(0.1, 0.6)), 3))

    valid = [r for r in reactivities if r != -999]
    print(f'Sequence length: {len(sequence)}; valid reactivities: {len(valid)}/{len(reactivities)}; mean {np.mean(valid):.3f}')
    print()
    unconstrained, constrained = compare_predictions(sequence, reactivities)
    print('\n=== Plot ===')
    plot_reactivity(sequence, constrained['mfe'], reactivities)
    print('\nDone.')
