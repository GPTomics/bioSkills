#!/usr/bin/env python3
'''
SHAPE-constrained RNA structure prediction.
Loads ShapeMapper2 reactivity profiles and uses them to constrain ViennaRNA folding.
'''
# Reference: ViennaRNA 2.6+, matplotlib 3.8+, numpy 1.26+, pandas 2.2+ | Verify API if version differs

import RNA
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def load_shape_profile(profile_file):
    '''Load ShapeMapper2 *_profile.txt output.'''
    df = pd.read_csv(profile_file, sep='\t')
    sequence = ''.join(df['Nucleotide'].tolist())
    reactivities = df['Reactivity_profile'].tolist()
    stderr = df['Std_err'].tolist() if 'Std_err' in df.columns else None
    return sequence, reactivities, stderr


def fold_unconstrained(sequence):
    '''Fold without any experimental constraints.'''
    fc = RNA.fold_compound(sequence)
    structure, mfe = fc.mfe()
    fc.pf()
    centroid, _ = fc.centroid()
    diversity = fc.mean_bp_distance()
    return {'mfe': structure, 'mfe_energy': mfe, 'centroid': centroid, 'diversity': diversity}


def fold_shape_constrained(sequence, reactivities, m=1.8, b=-0.6):
    '''
    Fold with SHAPE reactivity constraints.

    Deigan et al. (2009) pseudo-energy parameters:
    m=1.8, b=-0.6: standard for SHAPE (1M7, NAI, NMIA)
    m=1.1, b=-0.3: suggested for DMS-MaPseq
    '''
    shape_data = [-999.0] + [r if r != -999 else -999.0 for r in reactivities]

    fc = RNA.fold_compound(sequence)
    fc.sc_add_SHAPE_deigan(shape_data, m, b)
    structure, mfe = fc.mfe()

    fc2 = RNA.fold_compound(sequence)
    fc2.sc_add_SHAPE_deigan(shape_data, m, b)
    fc2.pf()
    centroid, _ = fc2.centroid()
    diversity = fc2.mean_bp_distance()

    return {'mfe': structure, 'mfe_energy': mfe, 'centroid': centroid, 'diversity': diversity}


def shape_agreement(structure, reactivities, low=0.3, high=0.7):
    '''
    Fraction of positions where SHAPE agrees with structure.

    Paired positions should have low reactivity (< 0.3).
    Unpaired positions should have high reactivity (> 0.7).
    Positions with moderate reactivity (0.3-0.7) or no data are excluded.
    '''
    agree, total = 0, 0
    for char, react in zip(structure, reactivities):
        if react == -999 or (react >= low and react <= high):
            continue
        total += 1
        paired = char in '()'
        if paired and react < low:
            agree += 1
        elif not paired and react > high:
            agree += 1
    return agree / total if total > 0 else 0.0


def compare_predictions(sequence, reactivities, output_prefix='shape_comparison'):
    '''Compare unconstrained and SHAPE-constrained predictions.'''
    unconstrained = fold_unconstrained(sequence)
    constrained = fold_shape_constrained(sequence, reactivities)

    bp_dist = RNA.bp_distance(unconstrained['mfe'], constrained['mfe'])

    unconstr_agreement = shape_agreement(unconstrained['mfe'], reactivities)
    constr_agreement = shape_agreement(constrained['mfe'], reactivities)

    print('=== Unconstrained ===')
    print(f'  MFE: {unconstrained["mfe"]} ({unconstrained["mfe_energy"]:.2f} kcal/mol)')
    print(f'  Centroid: {unconstrained["centroid"]}')
    print(f'  Ensemble diversity: {unconstrained["diversity"]:.2f}')
    print(f'  SHAPE agreement: {unconstr_agreement:.1%}')

    print('\n=== SHAPE-constrained ===')
    print(f'  MFE: {constrained["mfe"]} ({constrained["mfe_energy"]:.2f} kcal/mol)')
    print(f'  Centroid: {constrained["centroid"]}')
    print(f'  Ensemble diversity: {constrained["diversity"]:.2f}')
    print(f'  SHAPE agreement: {constr_agreement:.1%}')

    print(f'\nBP distance between predictions: {bp_dist}')

    return unconstrained, constrained


def plot_reactivity_with_structure(sequence, structure, reactivities, output_file='reactivity_structure.png'):
    '''Plot reactivity profile colored by predicted pairing status.'''
    n = len(sequence)
    positions = np.arange(1, n + 1)

    valid = np.array([r != -999 for r in reactivities])
    react_arr = np.array([r if r != -999 else 0 for r in reactivities])

    paired = np.array([c in '()' for c in structure])
    colors = np.where(paired, '#4169E1', '#FF4500')

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(max(8, n * 0.08), 6),
                                    gridspec_kw={'height_ratios': [3, 1]}, sharex=True)

    ax1.bar(positions[valid], react_arr[valid], color=colors[valid], width=1.0, edgecolor='none')
    ax1.axhline(y=0.3, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax1.axhline(y=0.7, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax1.set_ylabel('SHAPE reactivity')
    ax1.set_title('SHAPE reactivity (blue=paired, red=unpaired)')

    struct_binary = np.where(paired, 1.0, 0.0)
    ax2.fill_between(positions, struct_binary, step='mid', alpha=0.3, color='steelblue')
    ax2.set_ylabel('Paired')
    ax2.set_xlabel('Position')
    ax2.set_yticks([0, 1])
    ax2.set_yticklabels(['Unpaired', 'Paired'])
    ax2.set_xlim(0, n + 1)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f'Saved to {output_file}')


if __name__ == '__main__':
    profile_file = 'shapemapper_results/my_rna_profile.txt'

    try:
        sequence, reactivities, stderr = load_shape_profile(profile_file)
    except FileNotFoundError:
        print(f'Profile not found: {profile_file}')
        print('Using example data instead.')

        sequence = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'
        # Simulated SHAPE reactivities for tRNA
        np.random.seed(42)
        n = len(sequence)
        reactivities = []
        for i in range(n):
            if i < 7 or (i > 25 and i < 43) or i > 65:
                reactivities.append(round(np.random.uniform(0.0, 0.3), 3))
            elif (i >= 14 and i <= 21) or (i >= 32 and i <= 37):
                reactivities.append(round(np.random.uniform(0.5, 1.2), 3))
            else:
                reactivities.append(round(np.random.uniform(0.1, 0.6), 3))
        stderr = None

    print(f'Sequence length: {len(sequence)}')
    valid_react = [r for r in reactivities if r != -999]
    print(f'Valid reactivities: {len(valid_react)}/{len(reactivities)}')
    print(f'Mean reactivity: {np.mean(valid_react):.3f}')

    print('\n=== Comparing predictions ===')
    unconstrained, constrained = compare_predictions(sequence, reactivities)

    print('\n=== Generating plots ===')
    plot_reactivity_with_structure(sequence, constrained['mfe'], reactivities)

    print('\nDone.')
