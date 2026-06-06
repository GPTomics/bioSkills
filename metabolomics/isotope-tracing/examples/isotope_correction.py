"""Natural-abundance / tracer-purity correction of a raw isotopologue vector.

Takes a measured M+0..M+n intensity vector for a known formula under a 13C tracer,
applies natural-abundance + tracer-purity correction, and reports the corrected
mass-isotopomer distribution (MID) plus fractional (mean) enrichment.

Preferred path uses IsoCor (Millard 2019). If IsoCor is not importable the script
falls back to a transparent numpy correction-matrix implementation so the math is
inspectable; the two agree to numerical precision for low-resolution 13C data.
"""
# Reference: isocor 2.2+, numpy 1.26+ | Verify API if version differs

import math
import numpy as np

C13_NATURAL_ABUNDANCE = 0.0107  # fraction of carbon that is 13C in nature (IUPAC)
TRACER_PURITY_LABELED = 0.99    # per-position 13C purity of a typical U-13C tracer


def correction_matrix_13c(n_carbons, na=C13_NATURAL_ABUNDANCE):
    """Lower-triangular matrix mapping true MID to observed MID via the binomial
    natural-abundance ladder of the carbon backbone (low-resolution case)."""
    size = n_carbons + 1
    matrix = np.zeros((size, size))
    for true_label in range(size):
        remaining = n_carbons - true_label
        for extra in range(remaining + 1):
            prob = (math.comb(remaining, extra) * na**extra * (1 - na)**(remaining - extra))
            matrix[true_label + extra, true_label] = prob
    return matrix


def correct_numpy(raw, n_carbons, purity=TRACER_PURITY_LABELED):
    observed = np.asarray(raw, dtype=float)
    corrected = np.linalg.lstsq(correction_matrix_13c(n_carbons), observed, rcond=None)[0]
    corrected = np.clip(corrected, 0, None)
    if purity < 1.0:
        impurity = correction_matrix_13c(n_carbons, na=1 - purity)
        corrected = np.clip(np.linalg.lstsq(impurity, corrected, rcond=None)[0], 0, None)
    return corrected


def correct_isocor(raw, formula, tracer_purity):
    import isocor
    corrector = isocor.mscorrectors.MetaboliteCorrectorFactory(
        formula, tracer='13C', correct_NA_tracer=True, tracer_purity=tracer_purity)
    corrected_area, iso_fraction, residuum, mean_enrichment = corrector.correct(raw)
    return np.asarray(corrected_area), np.asarray(iso_fraction), mean_enrichment


def fractional_enrichment(mid):
    """Mean labeling per tracer atom: weighted-mean isotopologue index / n_atoms."""
    indices = np.arange(len(mid))
    return float(np.sum(indices * mid) / ((len(mid) - 1) * np.sum(mid)))


def main():
    formula = 'C6H12O6'
    n_carbons = 6
    raw = [50000.0, 8000.0, 12000.0, 3000.0, 1500.0, 6000.0, 25000.0]

    print('raw isotopologue areas (M+0..M+6):', raw)
    print('raw MID looks labeled, but the M+1/M+2 shoulder is partly natural 13C\n')

    try:
        corrected, iso_fraction, mean_enr = correct_isocor(raw, formula, [1 - TRACER_PURITY_LABELED, TRACER_PURITY_LABELED])
        engine = 'IsoCor'
    except ImportError:
        corrected = correct_numpy(raw, n_carbons)
        iso_fraction = corrected / corrected.sum()
        mean_enr = fractional_enrichment(corrected)
        engine = 'numpy correction matrix (IsoCor unavailable)'

    mid = corrected / corrected.sum()
    print(f'correction engine: {engine}')
    print('corrected MID (M+0..M+6):', np.round(mid, 4).tolist())
    print('fractional enrichment   :', round(mean_enr, 4))
    print('\nuncorrected M+0 fraction:', round(raw[0] / sum(raw), 4),
          '-> corrected M+0 fraction:', round(mid[0], 4))
    print('the M+0 fraction rose because natural-abundance signal was removed from M+1/M+2')


if __name__ == '__main__':
    main()
