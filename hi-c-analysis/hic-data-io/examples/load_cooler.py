'''Load a cooler at one resolution, fetch balanced pixels, and build a cooler from a matrix (vectorized).'''
# Reference: cooler 0.10+, bioframe 0.7+ | Verify API if version differs

import cooler
import bioframe
import numpy as np
import pandas as pd

BIN_SIZE = 10_000          # 10kb: a common Hi-C analysis resolution; the bin table defines pixel comparability
COARSEN_RES = [10_000, 50_000, 100_000, 500_000]   # integer-multiple ladder so HiGlass tiling stays valid


def load_and_inspect(mcool_path, resolution):
    '''Resolve to a single-resolution URI, report info, and fetch the balanced chr1 matrix.'''
    available = cooler.fileops.list_coolers(mcool_path)
    print(f'Available resolutions: {available}')
    clr = cooler.Cooler(f'{mcool_path}::/resolutions/{resolution}')   # never the bare .mcool
    print(f'Resolution: {clr.binsize} bp | chromosomes: {len(clr.chromnames)} | total contacts: {clr.info["sum"]:,}')

    balanced = 'weight' in clr.bins().columns
    print(f'Balanced (has weight column): {balanced}')
    if not balanced:
        print('Not balanced - matrix(balance=True) would raise; balance first (see matrix-operations).')
        return clr

    matrix = clr.matrix(balance=True).fetch('chr1')                  # cooler weight is multiplicative
    masked = np.isnan(matrix).all(axis=1).sum()                      # all-NaN rows = masked low-coverage bins, not a bug
    print(f'chr1 matrix shape: {matrix.shape} | masked bins: {masked} | max balanced: {np.nanmax(matrix):.4f}')
    return clr


def matrix_to_cooler(matrix, out_path, assembly='hg38'):
    '''Write an in-memory contact matrix to a cooler without a per-bin-pair loop.'''
    chromsizes = bioframe.fetch_chromsizes(assembly)                 # pin assembly + chrom naming up front
    bins = cooler.binnify(chromsizes, BIN_SIZE)

    upper = np.triu(matrix)                                          # cooler stores the upper triangle only
    i, j = np.nonzero(upper)                                         # vectorized; no O(n^2) Python loop
    pixels = pd.DataFrame({'bin1_id': i, 'bin2_id': j, 'count': upper[i, j]})
    cooler.create_cooler(out_path, bins, pixels, assembly=assembly, symmetric_upper=True)
    print(f'Wrote {len(pixels):,} pixels to {out_path}')


def coarsen(hires_path, out_mcool):
    '''Coarsen correctly: zoomify sums raw counts then re-runs ICE per resolution.'''
    cooler.zoomify_cooler(hires_path, out_mcool, resolutions=COARSEN_RES, chunksize=10_000_000)
    print(f'Zoomified to {COARSEN_RES} (raw-summed, re-balance per level)')


if __name__ == '__main__':
    clr = load_and_inspect('matrix.mcool', BIN_SIZE)
