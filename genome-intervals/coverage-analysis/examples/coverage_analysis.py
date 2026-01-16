#!/usr/bin/env python3
'''Demonstrate coverage analysis with pybedtools.'''

import pybedtools
import subprocess

# Create sample BED data for demonstration
regions_str = '''chr1\t100\t500\tregion1\t0\t+
chr1\t1000\t1500\tregion2\t0\t-
chr2\t200\t600\tregion3\t0\t+'''

reads_str = '''chr1\t150\t250\tread1\t0\t+
chr1\t200\t300\tread2\t0\t+
chr1\t250\t350\tread3\t0\t+
chr1\t300\t400\tread4\t0\t+
chr1\t1100\t1200\tread5\t0\t-
chr2\t250\t350\tread6\t0\t+
chr2\t300\t400\tread7\t0\t+'''

genome_str = '''chr1\t2000
chr2\t1000'''

regions = pybedtools.BedTool(regions_str, from_string=True)
reads = pybedtools.BedTool(reads_str, from_string=True)

# Write genome file
with open('test_genome.txt', 'w') as f:
    f.write(genome_str)

print('=== Input Data ===')
print(f'Regions: {regions.count()} intervals')
print(f'Reads: {reads.count()} intervals')

# Coverage: how much of each region is covered by reads
print('\n=== Coverage per Region ===')
coverage = regions.coverage(reads)
print('Format: chrom, start, end, name, score, strand, overlaps, bases_covered, region_length, fraction')
for interval in coverage:
    print(f'  {interval}')

# Extract just the important info
print('\n=== Coverage Summary ===')
for interval in regions.coverage(reads):
    fields = interval.fields
    region_name = fields[3]
    overlapping_reads = int(fields[6])
    bases_covered = int(fields[7])
    region_length = int(fields[8])
    fraction = float(fields[9])
    print(f'  {region_name}: {overlapping_reads} reads, {bases_covered}/{region_length} bases ({fraction:.1%} covered)')

# Coverage with counts only
print('\n=== Read Counts per Region ===')
counts = regions.coverage(reads, counts=True)
for interval in counts:
    fields = interval.fields
    print(f'  {fields[3]}: {fields[6]} overlapping reads')

# Generate coverage bedGraph from BED intervals
print('\n=== Coverage bedGraph (from reads) ===')
reads_sorted = reads.sort()
bedgraph = reads_sorted.genome_coverage(bg=True, g='test_genome.txt')
print('First few lines:')
for i, interval in enumerate(bedgraph):
    if i >= 5:
        print('  ...')
        break
    print(f'  {interval.chrom}\t{interval.start}\t{interval.end}\t{interval.fields[3]}')

# Mean coverage per region
print('\n=== Mean Coverage per Region ===')
mean_cov = regions.coverage(reads, mean=True)
for interval in mean_cov:
    fields = interval.fields
    print(f'  {fields[3]}: mean coverage = {float(fields[6]):.2f}x')

# Identify low-coverage regions
print('\n=== Low Coverage Regions (<50% covered) ===')
low_cov = regions.coverage(reads).filter(lambda x: float(x.fields[9]) < 0.5)
for interval in low_cov:
    print(f'  {interval.fields[3]}: {float(interval.fields[9]):.1%} covered')

# Cleanup
import os
os.remove('test_genome.txt')
pybedtools.cleanup()

print('\n=== Done ===')
