#!/usr/bin/env python3
'''Demonstrate core interval arithmetic operations with pybedtools.'''

import pybedtools

# Create sample BED files
peaks_str = '''chr1\t100\t200\tpeak1\t100\t+
chr1\t300\t400\tpeak2\t200\t+
chr1\t500\t600\tpeak3\t150\t+
chr2\t100\t200\tpeak4\t250\t-'''

genes_str = '''chr1\t150\t350\tgeneA\t0\t+
chr1\t550\t700\tgeneB\t0\t-
chr2\t50\t150\tgeneC\t0\t+'''

peaks = pybedtools.BedTool(peaks_str, from_string=True)
genes = pybedtools.BedTool(genes_str, from_string=True)

print('=== Original Data ===')
print(f'Peaks: {peaks.count()} intervals')
print(f'Genes: {genes.count()} intervals')

# Intersect - find overlapping regions
print('\n=== Intersect (overlapping portions) ===')
overlap = peaks.intersect(genes)
for i in overlap:
    print(f'  {i.chrom}:{i.start}-{i.end}')

# Intersect -u (report A entries that overlap B)
print('\n=== Intersect -u (peaks overlapping genes) ===')
peaks_in_genes = peaks.intersect(genes, u=True)
for i in peaks_in_genes:
    print(f'  {i.chrom}:{i.start}-{i.end} ({i.name})')

# Intersect -v (report A entries that DON'T overlap B)
print('\n=== Intersect -v (peaks NOT overlapping genes) ===')
peaks_outside = peaks.intersect(genes, v=True)
for i in peaks_outside:
    print(f'  {i.chrom}:{i.start}-{i.end} ({i.name})')

# Intersect -wa -wb (report both A and B fields)
print('\n=== Intersect -wa -wb (with gene info) ===')
with_info = peaks.intersect(genes, wa=True, wb=True)
for i in with_info:
    fields = i.fields
    print(f'  Peak: {fields[3]} overlaps Gene: {fields[9]}')

# Subtract
print('\n=== Subtract (remove gene regions from peaks) ===')
subtracted = peaks.subtract(genes)
for i in subtracted:
    print(f'  {i.chrom}:{i.start}-{i.end}')

# Merge
print('\n=== Merge (combine nearby intervals) ===')
# Create overlapping intervals for merge demo
merge_str = '''chr1\t100\t200
chr1\t150\t250
chr1\t240\t300
chr1\t500\t600'''
to_merge = pybedtools.BedTool(merge_str, from_string=True)
merged = to_merge.sort().merge()
print('Before merge:')
for i in to_merge:
    print(f'  {i.chrom}:{i.start}-{i.end}')
print('After merge:')
for i in merged:
    print(f'  {i.chrom}:{i.start}-{i.end}')

# Merge with distance tolerance
merged_d100 = to_merge.sort().merge(d=100)
print('After merge with d=100:')
for i in merged_d100:
    print(f'  {i.chrom}:{i.start}-{i.end}')

# Jaccard similarity
print('\n=== Jaccard Similarity ===')
jaccard = peaks.jaccard(genes)
print(f"  Intersection: {jaccard['intersection']} bp")
print(f"  Union: {jaccard['union']} bp")
print(f"  Jaccard index: {jaccard['jaccard']:.4f}")

# Cleanup temp files
pybedtools.cleanup()
print('\n=== Done ===')
