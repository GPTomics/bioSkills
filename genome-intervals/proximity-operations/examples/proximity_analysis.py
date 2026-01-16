#!/usr/bin/env python3
'''Demonstrate proximity operations with pybedtools.'''

import pybedtools

# Create sample data
peaks_str = '''chr1\t1000\t1100\tpeak1\t100\t+
chr1\t5000\t5200\tpeak2\t150\t+
chr1\t20000\t20500\tpeak3\t200\t-
chr2\t3000\t3100\tpeak4\t120\t+'''

genes_str = '''chr1\t500\t2000\tgeneA\t0\t+
chr1\t8000\t15000\tgeneB\t0\t+
chr1\t18000\t25000\tgeneC\t0\t-
chr2\t1000\t5000\tgeneD\t0\t+'''

# Create genome file for slop/flank
genome_str = '''chr1\t50000
chr2\t50000'''

peaks = pybedtools.BedTool(peaks_str, from_string=True)
genes = pybedtools.BedTool(genes_str, from_string=True)

# Write genome file
with open('test_genome.txt', 'w') as f:
    f.write(genome_str)

print('=== Input Data ===')
print('Peaks:')
for p in peaks:
    print(f'  {p.name}: {p.chrom}:{p.start}-{p.end}')
print('\nGenes:')
for g in genes:
    print(f'  {g.name}: {g.chrom}:{g.start}-{g.end}')

# Closest operation
print('\n=== Closest Gene to Each Peak ===')
closest = peaks.closest(genes, d=True)
for interval in closest:
    fields = interval.fields
    peak_name = fields[3]
    gene_name = fields[9]
    distance = fields[-1]
    print(f'  {peak_name} -> {gene_name} (distance: {distance})')

# Window operation
print('\n=== Genes Within 5kb of Peaks ===')
window_result = peaks.window(genes, w=5000)
for interval in window_result:
    fields = interval.fields
    peak_name = fields[3]
    gene_name = fields[9]
    print(f'  {peak_name} near {gene_name}')

# Slop operation (extend intervals)
print('\n=== Extended Peaks (+500bp each side) ===')
extended = peaks.slop(g='test_genome.txt', b=500)
for original, ext in zip(peaks, extended):
    print(f'  {original.name}: {original.start}-{original.end} -> {ext.start}-{ext.end}')

# Flank operation (get flanking regions, not original)
print('\n=== Upstream Flanks (500bp) ===')
flanks = peaks.flank(g='test_genome.txt', l=500, r=0, s=True)
for original, flank in zip(peaks, flanks):
    print(f'  {original.name} ({original.strand}): upstream flank {flank.start}-{flank.end}')

# Find peaks within distance of gene TSS
print('\n=== Create TSS and Find Nearby Peaks ===')
# Extract TSS from genes
tss_intervals = []
for gene in genes:
    if gene.strand == '+':
        tss_intervals.append((gene.chrom, gene.start, gene.start + 1, gene.name, 0, gene.strand))
    else:
        tss_intervals.append((gene.chrom, gene.end - 1, gene.end, gene.name, 0, gene.strand))
tss = pybedtools.BedTool(tss_intervals)

print('TSS positions:')
for t in tss:
    print(f'  {t.name}: {t.chrom}:{t.start}')

peaks_near_tss = peaks.window(tss, w=3000)
print('\nPeaks within 3kb of TSS:')
for interval in peaks_near_tss:
    fields = interval.fields
    print(f'  {fields[3]} near {fields[9]} TSS')

# Cleanup
import os
os.remove('test_genome.txt')
pybedtools.cleanup()

print('\n=== Done ===')
