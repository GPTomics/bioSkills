#!/usr/bin/env python3
'''Demonstrate pyBigWig for reading and creating bigWig files.'''

import pyBigWig
import numpy as np

# Create a sample bigWig file
print('=== Creating BigWig File ===')
bw = pyBigWig.open('sample.bw', 'w')

# Add chromosome sizes header
chroms = [('chr1', 10000), ('chr2', 8000)]
bw.addHeader(chroms)
print(f'Added header: {chroms}')

# Add entries for chr1
# Method 1: Variable width intervals
bw.addEntries(['chr1', 'chr1', 'chr1'],
              [0, 500, 1000],
              ends=[500, 1000, 2000],
              values=[1.5, 2.0, 3.5])
print('Added variable-width entries to chr1')

# Method 2: Fixed-width spans (more efficient for dense data)
bw.addEntries('chr1',
              3000,  # start position
              values=[4.0, 5.0, 6.0, 5.5, 4.5],
              span=100,
              step=100)
print('Added fixed-span entries to chr1 (3000-3500)')

# Add entries for chr2
bw.addEntries('chr2', [0, 1000, 2000],
              ends=[1000, 2000, 3000],
              values=[2.0, 3.0, 1.0])
print('Added entries to chr2')

bw.close()
print('BigWig file created: sample.bw')

# Read the bigWig file
print('\n=== Reading BigWig File ===')
bw = pyBigWig.open('sample.bw')

# File information
print(f'Is bigWig: {bw.isBigWig()}')
print(f'Chromosomes: {bw.chroms()}')
print(f'Header: {bw.header()}')

# Get values for a region
print('\n=== Get Values for Region ===')
values = bw.values('chr1', 0, 1000)
print(f'chr1:0-1000 values (first 10): {values[:10]}')
print(f'Mean: {np.nanmean(values):.2f}')
print(f'Max: {np.nanmax(values):.2f}')

# Get intervals (non-zero regions)
print('\n=== Get Intervals ===')
intervals = bw.intervals('chr1', 0, 4000)
print('chr1:0-4000 intervals:')
for start, end, val in intervals:
    print(f'  {start}-{end}: {val}')

# Statistics for regions
print('\n=== Statistics ===')
stats_types = ['mean', 'min', 'max', 'coverage', 'std', 'sum']
for stat in stats_types:
    result = bw.stats('chr1', 0, 2000, type=stat)
    print(f'chr1:0-2000 {stat}: {result[0]:.3f}' if result[0] else f'chr1:0-2000 {stat}: None')

# Binned statistics
print('\n=== Binned Statistics ===')
binned = bw.stats('chr1', 0, 4000, type='mean', nBins=8)
print(f'chr1:0-4000 in 8 bins (500bp each):')
for i, val in enumerate(binned):
    start = i * 500
    end = (i + 1) * 500
    print(f'  {start}-{end}: {val:.2f}' if val else f'  {start}-{end}: None')

bw.close()

# Cleanup
import os
os.remove('sample.bw')
print('\n=== Done (sample.bw cleaned up) ===')
