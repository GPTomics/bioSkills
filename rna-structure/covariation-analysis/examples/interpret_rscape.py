#!/usr/bin/env python3
'''
Turn R-scape covariation output into the three-way verdict a structure claim needs.

The mistake R-scape exists to prevent is reading "0 significant pairs" as "no structure".
A negative is only meaningful if the alignment had the POWER to detect covariation, so the
verdict depends on BOTH the number of significantly covarying pairs AND the statistical power.
The per-pair numbers (E-value <= target counts as significant; per-pair power) come from the
R-scape .cov / .power outputs; verify the exact column order with `head -20 <msa>.cov`.
'''
# Reference: R-scape 2.0+ | Verify API if version differs


def classify_structure_support(n_significant, expected_pairs, mean_power, evalue_target=0.05,
                               power_floor=0.10):
    '''
    Three-way verdict from R-scape numbers.

    n_significant: pairs with covariation E-value <= evalue_target
    expected_pairs: base pairs in the proposed/predicted structure
    mean_power: R-scape 'alignment power' = sum of per-pair power / nbpairs (mean per-pair power)
    power_floor: R-scape's own (explicitly arbitrary) 10% threshold for low- vs high-power (Rivas 2020)
    '''
    if n_significant > 0:
        return 'supports', f'{n_significant}/{expected_pairs} pairs covary significantly (E <= {evalue_target})'
    if mean_power >= power_floor:
        return 'rejects', f'adequate power (mean {mean_power:.2f}) but no significant covariation -> structure not supported'
    return 'cannot_infer', f'low power (mean {mean_power:.2f}); gather more diverse homologs before concluding'


def parse_cov_header(cov_file):
    '''Read the R-scape .cov header (nseq, alen, nbpairs) if present; tolerant of format drift.'''
    info = {}
    with open(cov_file) as f:
        for line in f:
            if not line.startswith('#'):
                break
            for key in ('nseq', 'alen', 'avgid', 'nbpairs'):
                if key in line:
                    after = line.split(key, 1)[1].split()
                    if after:
                        info[key] = after[0]
    return info


if __name__ == '__main__':
    # Representative scenarios (the numbers come from R-scape's .cov / .power summaries).
    scenarios = [
        ('Rfam SEED (deep, diverse)', dict(n_significant=18, expected_pairs=30, mean_power=0.55)),
        ('Proposed lncRNA structure', dict(n_significant=0, expected_pairs=45, mean_power=0.42)),
        ('Shallow alignment, 4 similar seqs', dict(n_significant=0, expected_pairs=25, mean_power=0.03)),
    ]
    for name, nums in scenarios:
        verdict, reason = classify_structure_support(**nums)
        print(f'{name}:')
        print(f'  verdict: {verdict}')
        print(f'  {reason}')
        print()
    print('Rule: only a POWERED negative ("rejects") argues against a structure; a low-power')
    print('negative ("cannot infer") says nothing -- the lncRNA case is "rejects", not "cannot infer".')
