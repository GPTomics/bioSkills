#!/usr/bin/env python3
'''Query and write a bigWig the right way: statistic = question, exact for results, NaN is not zero.'''
# Reference: pyBigWig 0.3.22+, numpy 1.26+ | Verify API if version differs

import pyBigWig
import numpy as np
import os

# Build a toy track: 1 Mb of flat background (1.0) with a narrow 200 bp summit (500) and a no-data gap.
# This is the canonical case where mean and max disagree and where gaps are NaN.
BACKGROUND = 1.0
SUMMIT = 500.0
SUMMIT_START, SUMMIT_END = 500_000, 500_200   # narrow tall feature
GAP_START, GAP_END = 800_000, 900_000         # 100 kb left as no-data (NaN, not 0)
CHROM_LEN = 1_000_000

bw = pyBigWig.open('demo.bw', 'w')
bw.addHeader([('chr1', CHROM_LEN)])   # ordered (name,length); MUST precede addEntries; default builds 10 zoom levels
starts, ends, vals = [], [], []
cursor = 0
for s, e, v in [(0, SUMMIT_START, BACKGROUND), (SUMMIT_START, SUMMIT_END, SUMMIT),
                (SUMMIT_END, GAP_START, BACKGROUND), (GAP_END, CHROM_LEN, BACKGROUND)]:
    starts.append(s); ends.append(e); vals.append(v)
bw.addEntries(['chr1'] * len(starts), starts, ends=ends, values=vals)   # entries in sorted (chrom,start) order
bw.close()

bw = pyBigWig.open('demo.bw')
print('isBigWig:', bw.isBigWig(), '| chroms:', bw.chroms(), '| zoom levels:', bw.header()['nLevels'])

# Trap 1: WHICH statistic. A single mean over the megabase buries the summit; max exposes it.
region = ('chr1', 0, CHROM_LEN)
print('\nstatistic = biological question (exact=True so the number is real, not a zoom approximation):')
print('  mean (default; dilutes the peak):', bw.stats(*region, type='mean', exact=True)[0])
print('  max  (peak height):              ', bw.stats(*region, type='max', exact=True)[0])
print('  sum  (total amount):             ', bw.stats(*region, type='sum', exact=True)[0])
print('  coverage (fraction with data):   ', bw.stats(*region, type='coverage', exact=True)[0])

# Trap 2: granularity. One bin hides the summit; a max profile keeps it visible.
profile = bw.stats(*region, type='max', nBins=10)
print('\n10-bin max profile (the summit-containing bin stands out):')
print('  ', [round(x, 1) if x is not None else None for x in profile])

# Trap 3: NaN is not zero. Three gap-handling choices, three different numbers.
v = bw.values(*region, numpy=True)   # list by default; numpy=True -> ndarray with nan for the gap
print('\ngap handling over a track that is partly no-data:')
print('  np.mean      (poisons to NaN):  ', np.mean(v))
print('  np.nanmean   (covered-only=mean):', round(float(np.nanmean(v)), 4))
print('  nan_to_num   (gaps as zero=mean0):', round(float(np.nan_to_num(v).mean()), 4))

# Raw stored runs (unresampled truth for a region)
print('\nintervals() around the summit:', bw.intervals('chr1', SUMMIT_START - 100, SUMMIT_END + 100))

bw.close()
os.remove('demo.bw')
print('\nDone (demo.bw cleaned up)')
