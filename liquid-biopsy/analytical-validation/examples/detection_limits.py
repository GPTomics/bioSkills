'''Liquid-biopsy detection-limit math: genome equivalents, Poisson detection, and an LoD95 probit fit.

Treats a ctDNA assay as a molecule-counting experiment at the Poisson edge: input mass sets the
sampling ceiling, error suppression sets the VAF floor, and the achieved LoD is the worse of the two.
'''
# Reference: numpy 1.26+, scipy 1.12+, statsmodels 0.14+ | Verify API if version differs

import numpy as np
import statsmodels.api as sm
from scipy.stats import poisson, binom, norm

GE_PER_NG = 330  # haploid ~3.3 pg -> strict 1 ng / 3.3 pg = 303; 330 is the diploid-6.6 pg/rounding convention


def genome_equivalents(input_ng):
    return input_ng * GE_PER_NG


def detection_probability(input_ng, vaf, min_mutant_molecules=1):
    '''P(>= min_mutant_molecules present) under Poisson(lambda = GE * VAF).'''
    lam = genome_equivalents(input_ng) * vaf
    return float(poisson.sf(min_mutant_molecules - 1, lam))


def ge_for_sampling_detection(vaf, target_lambda=3.0):
    '''GE so lambda >= 3 -> ~95% chance the mutant molecule is present at all (1 - e^-3 = 0.95).'''
    return target_lambda / vaf


def limit_of_blank(blank_signals):
    '''CLSI EP17 LoB = mean + 1.645*SD; one-sided 95th percentile of analyte-free blanks.'''
    blank_signals = np.asarray(blank_signals, dtype=float)
    return blank_signals.mean() + 1.645 * blank_signals.std(ddof=1)


def lod95_probit(vaf_levels, detected):
    '''Probit fit of detection (0/1) on log10(VAF) - CLSI EP17 fits on log concentration,
    which respects the saturating detection curve; returns the VAF where P(detect) = 0.95.
    The dilution series must bracket the 0.95 crossing: all-detected upper levels cause
    near-complete separation and an unstable slope.'''
    log_vaf = np.log10(np.asarray(vaf_levels, dtype=float))
    y = np.asarray(detected, dtype=float)
    X = sm.add_constant(log_vaf)
    fit = sm.GLM(y, X, family=sm.families.Binomial(link=sm.families.links.Probit())).fit()
    intercept, slope = fit.params
    return float(10 ** ((norm.ppf(0.95) - intercept) / slope))


def panel_detection_probability(input_ng, vaf, n_loci, min_loci_positive=2):
    '''P(>= min_loci_positive of n_loci detected); >=2-of-N is the Signatera-style positivity rule.'''
    per_locus = detection_probability(input_ng, vaf)
    return float(binom.sf(min_loci_positive - 1, n_loci, per_locus))


def simulate_dilution_series(true_lod_vaf, input_ng, levels, replicates, seed=0):
    '''Simulate binary detection across a contrived dilution series for a probit LoD95 fit.'''
    rng = np.random.default_rng(seed)
    vaf_out, det_out = [], []
    for vaf in levels:
        lam = genome_equivalents(input_ng) * vaf
        present = rng.poisson(lam, size=replicates) >= 1
        vaf_out.extend([vaf] * replicates)
        det_out.extend(present.astype(int).tolist())
    return np.asarray(vaf_out), np.asarray(det_out)


def main():
    print('Poisson sampling ceiling')
    for ng in (3.3, 10, 30):
        print(f'  {ng:>5} ng = {genome_equivalents(ng):>6.0f} GE | P(detect 0.1% VAF) = {detection_probability(ng, 0.001):.3f}')
    print(f'  GE for 95% sampling-detection of a 0.01% variant: {ge_for_sampling_detection(1e-4):.0f} (~{ge_for_sampling_detection(1e-4)/GE_PER_NG:.0f} ng)')

    print('\nPer-locus vs panel-integrated detection at 30 ng, 1e-4 VAF')
    print(f'  single locus      P(detect) = {detection_probability(30, 1e-4):.3f}')
    for n in (16, 48):
        print(f'  {n:>2}-variant panel  P(>=2 detected) = {panel_detection_probability(30, 1e-4, n):.3f}')

    print('\nLoD95 from a simulated dilution series (true sampling LoD ~ where lambda=3)')
    levels = [5e-4, 1e-3, 2e-3, 5e-3, 1e-2]
    vaf, det = simulate_dilution_series(true_lod_vaf=1e-3, input_ng=10, levels=levels, replicates=24)
    blanks = np.random.default_rng(1).normal(0.0002, 0.00005, size=20)  # blank background-error signal
    print(f'  LoB (blank 95th pct VAF) = {limit_of_blank(blanks):.5f}')
    print(f'  LoD95 (probit, detection vs VAF) = {lod95_probit(vaf, det):.5f}')


if __name__ == '__main__':
    main()
