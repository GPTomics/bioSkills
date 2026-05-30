#!/usr/bin/env python3
'''Core interval arithmetic with pybedtools, plus the engine-equivalence note for pyranges/bioframe.'''
# Reference: pybedtools 0.10+, bedtools 2.31+ | Verify API if version differs

import pybedtools

MERGE_DIST = 100   # merge intervals within 100 bp into one (replicate-merge convention)

peaks_str = '''chr1\t100\t200\tpeak1\t100\t+
chr1\t300\t400\tpeak2\t200\t+
chr1\t500\t600\tpeak3\t150\t+
chr2\t100\t200\tpeak4\t250\t-'''

genes_str = '''chr1\t150\t350\tgeneA\t0\t+
chr1\t550\t700\tgeneB\t0\t-
chr2\t50\t150\tgeneC\t0\t+'''

peaks = pybedtools.BedTool(peaks_str, from_string=True)
genes = pybedtools.BedTool(genes_str, from_string=True)

print('peaks:', peaks.count(), 'genes:', genes.count())

print('\n=== intersect -u (whole A once if it overlaps any B) ===')
for i in peaks.intersect(genes, u=True):
    print(' ', i.name, f'{i.chrom}:{i.start}-{i.end}')

print('\n=== intersect -v (A with no overlap) ===')
for i in peaks.intersect(genes, v=True):
    print(' ', i.name)

print('\n=== intersect -wa -wb (join) ===')
for i in peaks.intersect(genes, wa=True, wb=True):
    print(' ', i.fields[3], 'overlaps', i.fields[9])

print('\n=== subtract (clip) vs subtract -A (drop whole A) ===')
print('  clipped:', [f'{i.chrom}:{i.start}-{i.end}' for i in peaks.subtract(genes)])
print('  dropped-any-overlap:', [i.name for i in peaks.subtract(genes, A=True)])

print('\n=== merge requires prior sort ===')
for i in peaks.sort().merge(d=MERGE_DIST, c='4,5', o='distinct,sum'):
    print(' ', i.fields)

print('\n=== map: aggregate a B column onto each sorted A interval ===')
for i in genes.sort().map(peaks.sort(), c=5, o='mean'):
    print(' ', i.name, '->', i.fields[-1])

print('\n=== jaccard is a similarity scalar, NOT a significance test ===')
j = peaks.jaccard(genes)
print(f"  jaccard={j['jaccard']:.4f} intersection={j['intersection']}bp union={j['union']}bp")

# Engine equivalence: pyranges and bioframe compute IDENTICAL geometry on the same
# 0-based half-open intervals. Porting bugs are about default strand handling and
# return shape (pyranges overlap vs join vs intersect; bioframe overlap how=), not
# the arithmetic. Check pyranges.__version__ for the v0/v1 API split before porting.

pybedtools.cleanup()
print('\n=== done ===')
