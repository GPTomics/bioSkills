#!/usr/bin/env python3
'''Genome-wide selection scan: Hudson FST, windowed Tajima's D and pi, iHS, and XP-EHH.

Demonstrates the load-bearing decisions of a sweep scan: empirical genome-wide outliers
instead of absolute cutoffs, ratio-of-averages FST after an MAF filter, the standardization
split (iHS within derived-allele-frequency bins, XP-EHH genome-wide), and intersecting
orthogonal signals. All outputs go to a caller-supplied directory (default a fresh temp dir)
so nothing is written to the current directory.
Usage: ./selection_scan.py <data.vcf.gz> [pop1.txt] [pop2.txt] [output_dir]
'''
# Reference: scikit-allel 1.3+, numpy 1.26+ | Verify API if version differs

import os
import sys
import tempfile
import allel
import numpy as np

WINDOW_SIZE = 100000          # bp; scale to the species' LD decay (smaller for high-recombination genomes)
WINDOW_STEP = 50000           # 50% overlap so adjacent windows share signal
MAF_FLOOR = 0.05              # rare variants deflate FST (Bhatia 2013) and make EHH undefined on singletons
OUTLIER_PCT = 99             # top 1% empirical tail; the bulk absorbs the shared demographic history
SCORE_FLAG = 2.0             # |standardized iHS/XP-EHH| > 2 ~ top 2.5% of N(0,1); a convention, not a calibrated p


def load_pop_index(path, samples):
    with open(path) as fh:
        wanted = {line.strip() for line in fh if line.strip()}
    return [i for i, s in enumerate(samples) if s in wanted]


def windowed_fst(pos, ac1, ac2, start, stop):
    '''Hudson FST: MAF-filtered, mean reported as ratio-of-averages not mean of per-SNP ratios.'''
    maf = np.minimum(ac1.to_frequencies()[:, 1], ac2.to_frequencies()[:, 1])
    keep = (maf > MAF_FLOOR) & (1 - maf > MAF_FLOOR)
    num, den = allel.hudson_fst(ac1[keep], ac2[keep])
    fst_mean = np.nansum(num) / np.nansum(den)
    # Pass an explicit start/stop so the window grid matches the Tajima's D grid for index-aligned intersection.
    fst_win, windows, _ = allel.windowed_hudson_fst(pos[keep], ac1[keep], ac2[keep],
                                                    size=WINDOW_SIZE, start=start, stop=stop, step=WINDOW_STEP)
    return fst_mean, fst_win, windows


def standardized_ihs(h, pos):
    '''iHS standardized WITHIN derived-allele-count bins - the raw score is frequency-dependent.'''
    ac = h.count_alleles()
    flt = (ac[:, 0] > 1) & (ac[:, 1] > 1)
    h_flt, pos_flt, ac_flt = h.compress(flt, axis=0), pos[flt], ac.compress(flt, axis=0)
    raw = allel.ihs(h_flt, pos_flt, min_maf=MAF_FLOOR, include_edges=True)
    std, _ = allel.standardize_by_allele_count(raw, ac_flt[:, 1])
    return pos_flt, std


def standardized_xpehh(h1, h2, pos):
    '''XP-EHH standardized with a plain GENOME-WIDE z-score - NOT the iHS frequency-bin rule.'''
    raw = allel.xpehh(h1, h2, pos, include_edges=True)
    return allel.standardize(raw)


def main():
    if len(sys.argv) < 2:
        sys.exit('Usage: selection_scan.py <data.vcf.gz> [pop1.txt] [pop2.txt] [output_dir]')

    vcf_file = sys.argv[1]
    pop1_file = sys.argv[2] if len(sys.argv) > 2 else None
    pop2_file = sys.argv[3] if len(sys.argv) > 3 else None
    outdir = sys.argv[4] if len(sys.argv) > 4 else tempfile.mkdtemp(prefix='selscan_')
    os.makedirs(outdir, exist_ok=True)
    print(f'Outputs -> {outdir}')

    callset = allel.read_vcf(vcf_file)
    gt = allel.GenotypeArray(callset['calldata/GT'])
    pos = callset['variants/POS']
    samples = list(callset['samples'])

    ac = gt.count_alleles()
    flt = ac.is_segregating() & (ac.max_allele() == 1)
    gt, pos, ac = gt.compress(flt, axis=0), pos[flt], ac.compress(flt, axis=0)
    print(f'Segregating biallelic SNPs: {gt.n_variants}')

    # Fix a common genomic window grid so FST and Tajima windows align by coordinate, not by chance.
    region_start, region_stop = int(pos.min()), int(pos.max())
    tajd, tajd_win, _ = allel.windowed_tajima_d(pos, ac, size=WINDOW_SIZE, start=region_start, stop=region_stop, step=WINDOW_STEP)
    pi, pi_win, _, _ = allel.windowed_diversity(pos, ac, size=WINDOW_SIZE, start=region_start, stop=region_stop, step=WINDOW_STEP)
    print(f"Windowed Tajima's D and pi over {len(tajd)} windows (interpret vs the genome-wide distribution)")

    if pop1_file and pop2_file:
        pop1_idx = load_pop_index(pop1_file, samples)
        pop2_idx = load_pop_index(pop2_file, samples)
        ac_sub = gt.count_alleles_subpops({'p1': pop1_idx, 'p2': pop2_idx})
        fst_mean, fst_win, fst_windows = windowed_fst(pos, ac_sub['p1'], ac_sub['p2'], region_start, region_stop)
        print(f'Mean Hudson FST (ratio-of-averages): {fst_mean:.4f}')

        # Same window grid, so window i is the same genomic interval in both scans.
        fst_tail = fst_win > np.nanpercentile(fst_win, OUTLIER_PCT)
        depressed = tajd < np.nanpercentile(tajd, 100 - OUTLIER_PCT)
        credible = np.where(fst_tail & depressed)[0]
        print(f'Windows in BOTH the FST tail and the low-Tajima tail (orthogonal signals): {len(credible)}')

        h = gt.to_haplotypes()
        pop1_hap = np.concatenate([[2 * i, 2 * i + 1] for i in pop1_idx]) if pop1_idx else np.array([], dtype=int)
        pop2_hap = np.concatenate([[2 * i, 2 * i + 1] for i in pop2_idx]) if pop2_idx else np.array([], dtype=int)
        xpehh_std = standardized_xpehh(h.take(pop1_hap, axis=1), h.take(pop2_hap, axis=1), pos)
        print(f'XP-EHH (completed sweeps) |z|>{SCORE_FLAG}: {(np.abs(xpehh_std) > SCORE_FLAG).sum()} SNPs')

        np.savetxt(os.path.join(outdir, 'fst_windows.tsv'),
                   np.column_stack([fst_windows[:, 0], fst_windows[:, 1], fst_win]),
                   header='start\tend\tfst', delimiter='\t', comments='')

    h = gt.to_haplotypes()
    ihs_pos, ihs_std = standardized_ihs(h, pos)
    print(f'iHS (incomplete sweeps) |std|>{SCORE_FLAG}: {(np.abs(ihs_std) > SCORE_FLAG).sum()} SNPs '
          f'(a null iHS does NOT mean no selection - completed sweeps need XP-EHH)')
    np.savetxt(os.path.join(outdir, 'ihs_scores.tsv'),
               np.column_stack([ihs_pos, ihs_std]), header='pos\tihs_std', delimiter='\t', comments='')

    print(f'Done. Tables written under {outdir}')


if __name__ == '__main__':
    main()
