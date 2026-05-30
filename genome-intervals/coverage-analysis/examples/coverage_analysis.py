#!/usr/bin/env python3
'''Interpret a coverage distribution: median, breadth curve, and evenness from a mosdepth dist file.'''
# Reference: mosdepth 0.3+, numpy 1.26+ | Verify API if version differs
# The point of this script: stop reporting a number, start reporting a curve.
# mosdepth global.dist.txt rows are: chrom, depth, proportion_of_bases_at_least_this_depth
# (the 'total' chrom is the whole-genome cumulative curve).

import numpy as np

BREADTH_THRESHOLDS = [1, 10, 20, 30]   # depths to report breadth at; 20x ~ germline SNV callability floor
SKEW_FLAG = 1.2                        # mean/median above this signals a right-skewed (tailed) distribution


def read_global_dist(path, chrom='total'):
    rows = [line.split() for line in open(path) if line.split()[0] == chrom]
    depths = np.array([int(r[1]) for r in rows])
    cum_ge = np.array([float(r[2]) for r in rows])
    order = np.argsort(depths)
    return depths[order], cum_ge[order]


def median_from_cumulative(depths, cum_ge):
    below_half = depths[cum_ge >= 0.5]
    return int(below_half.max()) if below_half.size else 0


def mean_from_cumulative(depths, cum_ge):
    per_depth = -np.diff(np.append(cum_ge, 0.0))   # fraction of bases exactly at each depth
    return float(np.sum(depths * per_depth))


def breadth_at(depths, cum_ge, threshold):
    at = cum_ge[depths == threshold]
    return float(at[0]) if at.size else float(cum_ge[depths >= threshold].max() if (depths >= threshold).any() else 0.0)


def evenness_cv(depths, cum_ge):
    per_depth = -np.diff(np.append(cum_ge, 0.0))
    mean = np.sum(depths * per_depth)
    var = np.sum(((depths - mean) ** 2) * per_depth)
    return float(np.sqrt(var) / mean) if mean else float('nan')


def report(path):
    depths, cum_ge = read_global_dist(path)
    median = median_from_cumulative(depths, cum_ge)
    mean = mean_from_cumulative(depths, cum_ge)
    cv = evenness_cv(depths, cum_ge)
    print(f'median depth: {median}x   mean depth: {mean:.1f}x   CV: {cv:.2f}')
    if median and mean / median > SKEW_FLAG:
        print(f'WARNING: mean/median = {mean/median:.2f} > {SKEW_FLAG} -- right-skewed; mean overstates typical depth (check dups/repeats/rDNA)')
    print('breadth curve:')
    for t in BREADTH_THRESHOLDS:
        print(f'  >= {t:2d}x: {breadth_at(depths, cum_ge, t):.1%}')


def demo():
    '''Synthetic 'total' dist: even ~30x library with a small right tail.'''
    lines = []
    rng = np.random.default_rng(0)
    sample = np.clip(rng.poisson(30, 200000), 0, None)
    sample = np.append(sample, rng.integers(800, 1200, 400))   # rDNA-style right tail inflating the mean
    maxd = int(sample.max())
    n = sample.size
    for d in range(maxd + 1):
        prop = float(np.mean(sample >= d))
        lines.append(f'total\t{d}\t{prop:.6f}\n')
    open('demo.mosdepth.global.dist.txt', 'w').writelines(lines)
    return 'demo.mosdepth.global.dist.txt'


if __name__ == '__main__':
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else demo()
    print(f'=== {path} ===')
    report(path)
