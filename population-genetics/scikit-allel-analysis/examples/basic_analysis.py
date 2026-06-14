#!/usr/bin/env python3
'''Decision-grade scikit-allel analysis: per-base diversity WITH an accessibility mask and
genome-wide FST as a ratio-of-sums (not a mean of per-SNP FST). Reads a VCF (or builds a small
synthetic callset when none is given) and writes a single TSV summary to a caller-supplied dir
(default: a fresh temp dir) so nothing is written to the current directory.'''
# Reference: scikit-allel 1.3.13+, numpy 1.26+ | Verify API if version differs

import os
import sys
import csv
import tempfile
import numpy as np
import allel


def load_or_simulate(vcf_file):
    '''Return (gt, pos, is_accessible, pop1_idx, pop2_idx). Simulates if no VCF is supplied.'''
    if vcf_file and os.path.exists(vcf_file):
        callset = allel.read_vcf(vcf_file, fields=['samples', 'calldata/GT', 'variants/POS'])
        gt = allel.GenotypeArray(callset['calldata/GT'])
        pos = callset['variants/POS']
        n = gt.n_samples
        half = n // 2
        return gt, pos, None, np.arange(half), np.arange(half, n)

    rng = np.random.default_rng(0)
    n_variants, n_per_pop, span = 600, 20, 1_000_000
    # Callable genome: two accessible flanks (300kb total) with a large inaccessible block in the
    # middle, so accessible bp is ~30% of the span. Variants sit only where callable, so the masked
    # numerator equals the naive one and the only difference is the denominator -> the deflation is real.
    is_accessible = np.zeros(span + 1, dtype=bool)
    is_accessible[1:150_001] = True
    is_accessible[850_000:span + 1] = True
    callable_pos = np.flatnonzero(is_accessible)
    pos = np.sort(rng.choice(callable_pos, size=n_variants, replace=False)).astype('i4')
    # Two populations with a frequency shift so FST is non-trivial.
    af1 = rng.uniform(0.05, 0.95, n_variants)
    af2 = np.clip(af1 + rng.normal(0, 0.15, n_variants), 0.01, 0.99)
    g1 = rng.binomial(2, af1[:, None], size=(n_variants, n_per_pop))
    g2 = rng.binomial(2, af2[:, None], size=(n_variants, n_per_pop))
    n012 = np.concatenate([g1, g2], axis=1)
    # Expand 012 dosages into a (n_variants, n_samples, 2) GenotypeArray.
    gt_data = np.stack([(n012 > 0).astype('i1'), (n012 > 1).astype('i1')], axis=-1)
    gt = allel.GenotypeArray(gt_data)
    return gt, pos, is_accessible, np.arange(n_per_pop), np.arange(n_per_pop, 2 * n_per_pop)


def main():
    vcf_file = sys.argv[1] if len(sys.argv) > 1 else None
    outdir = sys.argv[2] if len(sys.argv) > 2 else tempfile.mkdtemp()
    os.makedirs(outdir, exist_ok=True)

    gt, pos, is_accessible, pop1_idx, pop2_idx = load_or_simulate(vcf_file)

    ac = gt.count_alleles()
    flt = ac.is_segregating() & (ac.max_allele() == 1)
    gt, pos = gt.compress(flt, axis=0), pos[flt]
    ac = gt.count_alleles()

    subpops = {'pop1': list(pop1_idx), 'pop2': list(pop2_idx)}
    ac_subpops = gt.count_alleles_subpops(subpops)
    ac1, ac2 = ac_subpops['pop1'], ac_subpops['pop2']

    # Span the FULL region (including the inaccessible block) so the naive denominator is the span,
    # not just the variant range; otherwise the deflation the mask corrects for is invisible.
    if is_accessible is not None:
        start, stop = 1, len(is_accessible) - 1
    else:
        start, stop = None, None

    # Per-base diversity WITH the accessibility mask (the #1 silent bug is omitting it).
    pi_masked = allel.sequence_diversity(pos, ac, start=start, stop=stop, is_accessible=is_accessible)
    pi_naive = allel.sequence_diversity(pos, ac, start=start, stop=stop)   # divides by total span -> deflated
    theta_w = allel.watterson_theta(pos, ac, start=start, stop=stop, is_accessible=is_accessible)
    taj_d = allel.tajima_d(ac, pos=pos)

    # FST as a ratio-of-sums, NOT mean(per_snp_fst). Hudson is robust to unequal n (Bhatia 2013).
    fst_hudson, se, vb, vj = allel.average_hudson_fst(ac1, ac2, blen=50)   # blen > LD scale in real data
    num, den = allel.hudson_fst(ac1, ac2)
    fst_ratio = float(np.sum(num) / np.sum(den))
    per_snp = num / den
    fst_mean_wrong = float(np.nanmean(per_snp))   # the biased mean-of-ratios, for contrast only

    summary = [
        ('n_variants', gt.n_variants),
        ('n_samples', gt.n_samples),
        ('pi_with_accessibility_mask', round(pi_masked, 8)),
        ('pi_without_mask_deflated', round(pi_naive, 8)),
        ('watterson_theta', round(theta_w, 8)),
        ('tajima_d', round(float(taj_d), 5)),
        ('fst_hudson_ratio_of_sums', round(fst_ratio, 5)),
        ('fst_hudson_jackknife_se', round(float(se), 5)),
        ('fst_mean_of_per_snp_BIASED', round(fst_mean_wrong, 5)),
    ]

    out_path = os.path.join(outdir, 'allel_summary.tsv')
    with open(out_path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['statistic', 'value'])
        writer.writerows(summary)

    print(f'Outputs -> {outdir}')
    for k, v in summary:
        print(f'{k}\t{v}')
    print(f'Wrote {out_path}')


if __name__ == '__main__':
    main()
