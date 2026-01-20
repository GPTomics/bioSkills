'''Balance a Hi-C matrix'''

import cooler
import cooltools
import numpy as np

clr = cooler.Cooler('matrix.cool')

print('Before balancing:')
print(f'  Has weights: {"weight" in clr.bins().columns}')

print('\nBalancing (ICE)...')
cooler.balance_cooler('matrix.cool', store=True, cis_only=True)

clr = cooler.Cooler('matrix.cool')
weights = clr.bins()['weight'][:]
print(f'\nAfter balancing:')
print(f'  Weight range: {weights.min():.4f} - {weights.max():.4f}')
print(f'  NaN weights: {np.isnan(weights).sum()}')

print('\nComputing expected values...')
expected = cooltools.expected_cis(clr, ignore_diags=2)
print(expected.head(10))

balanced = clr.matrix(balance=True).fetch('chr1')
print(f'\nBalanced chr1 matrix:')
print(f'  Shape: {balanced.shape}')
print(f'  Range: {np.nanmin(balanced):.6f} - {np.nanmax(balanced):.6f}')
