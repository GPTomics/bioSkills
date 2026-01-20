'''Load and explore a Hi-C cooler file'''

import cooler
import numpy as np

clr = cooler.Cooler('matrix.mcool::resolutions/10000')

print('=== Hi-C Matrix Info ===')
print(f'Resolution: {clr.binsize} bp')
print(f'Chromosomes: {len(clr.chromnames)}')
print(f'Total bins: {clr.info["nbins"]}')
print(f'Total contacts: {clr.info["sum"]:,}')

bins = clr.bins()[:]
print(f'\nBin table shape: {bins.shape}')
print(bins.head())

print('\nExtracting chr1 matrix...')
matrix = clr.matrix(balance=True).fetch('chr1')
print(f'Matrix shape: {matrix.shape}')
print(f'Non-zero entries: {np.count_nonzero(~np.isnan(matrix)):,}')
print(f'Max value: {np.nanmax(matrix):.4f}')
