'''Target-decoy FDR and q-values from a concatenated-search PSM table.

Demonstrates the load-bearing idea: a raw PSM score is meaningless in isolation;
the actionable number is a list-level q-value (or a per-PSM PEP). This uses the
CONCATENATED target-decoy competition estimator FDR = decoys/targets. Separate
target/decoy searches would instead need the Elias-Gygi 2x-decoy form or the
distinct mix-max estimator (Keich 2015) -- do not mix them up.
Self-contained: builds a synthetic table so no input files are needed.'''
# Reference: numpy 1.26+, pandas 2.2+ | Verify API if version differs
import numpy as np
import pandas as pd

DECOY_PREFIXES = ('DECOY_', 'REV_', 'XXX_')
TARGET_FDR = 0.01   # 1% list-level FDR, the community standard for peptide IDs


def build_demo_table(n_target=2000, n_decoy=2000, seed=0):
    '''True targets score high, random targets and decoys score low and overlap.'''
    rng = np.random.default_rng(seed)
    true_targets = rng.normal(3.5, 0.8, n_target // 2)
    random_targets = rng.normal(1.0, 0.8, n_target - n_target // 2)
    decoys = rng.normal(1.0, 0.8, n_decoy)
    proteins = ['TARGET'] * n_target + [f'{DECOY_PREFIXES[0]}prot'] * n_decoy
    scores = np.concatenate([true_targets, random_targets, decoys])
    return pd.DataFrame({'protein': proteins, 'score': scores})


def add_qvalues(psms):
    '''Concatenated competition: rank by score, FDR = cumulative decoys / targets,
    then take the running minimum from the bottom to make q-values monotone.'''
    psms = psms.copy()
    psms['is_decoy'] = psms['protein'].str.startswith(DECOY_PREFIXES)
    psms = psms.sort_values('score', ascending=False).reset_index(drop=True)
    targets = (~psms['is_decoy']).cumsum()
    decoys = psms['is_decoy'].cumsum()
    psms['fdr'] = decoys / targets.clip(lower=1)
    psms['qvalue'] = psms['fdr'][::-1].cummin()[::-1]
    return psms


def add_pep(psms, bandwidth=0.3):
    '''PEP (local FDR) is per-PSM: the decoy density over total density at a score.
    Estimated here by a simple kernel-smoothed decoy fraction in a score window.'''
    psms = psms.copy()
    s = psms['score'].to_numpy()
    is_decoy = psms['is_decoy'].to_numpy().astype(float)
    pep = np.empty(len(s))
    for i, score in enumerate(s):
        w = np.exp(-((s - score) ** 2) / (2 * bandwidth ** 2))
        decoy_density = (is_decoy * w).sum()
        target_density = ((1 - is_decoy) * w).sum()
        pep[i] = min(1.0, decoy_density / max(target_density, 1e-9))   # local false-target rate, concatenated competition (no x2)
    psms['pep'] = pep
    return psms


def main():
    psms = add_qvalues(build_demo_table())
    psms = add_pep(psms)
    n_t = (~psms['is_decoy']).sum()
    n_d = psms['is_decoy'].sum()
    print(f'Targets: {n_t}, Decoys: {n_d} (concatenated 1:1 search)')

    kept = psms[(psms['qvalue'] <= TARGET_FDR) & (~psms['is_decoy'])]
    print(f'Target PSMs at q <= {TARGET_FDR}: {len(kept)}')
    print(f'Worst PEP inside the {TARGET_FDR:.0%}-FDR list: {kept["pep"].max():.3f}')

    strict = psms[(psms['pep'] <= TARGET_FDR) & (~psms['is_decoy'])]
    print(f'Target PSMs at PEP <= {TARGET_FDR}: {len(strict)} '
          f'(per-PSM cutoff is far stricter than the same q-value)')


if __name__ == '__main__':
    main()
