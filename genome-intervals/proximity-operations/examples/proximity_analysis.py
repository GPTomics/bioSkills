'''Proximity operations with pybedtools, done honestly.

Geometry (nearest TSS) is not biology (which gene an element regulates). For distal
elements nearest-gene is wrong the majority of the time (Fulco 2019); this script
computes the geometry correctly and flags distal hits for contact/QTL validation.
'''
# Reference: pybedtools 0.10+, bedtools 2.31+ | Verify API if version differs

import pybedtools

DISTAL_FLAG = 50000   # bp; |dist| beyond this is beyond promoter scale -> route to ABC/PCHi-C/eQTL, do not trust nearest

peaks_str = 'chr1\t1000\t1100\tpeak1\t100\t+\nchr1\t5000\t5200\tpeak2\t150\t+\nchr1\t20000\t20500\tpeak3\t200\t-\nchr2\t3000\t3100\tpeak4\t120\t+'
genes_str = 'chr1\t500\t2000\tgeneA\t0\t+\nchr1\t8000\t15000\tgeneB\t0\t+\nchr1\t18000\t25000\tgeneC\t0\t-\nchr2\t1000\t5000\tgeneD\t0\t+'
genome_str = 'chr1\t50000\nchr2\t50000'

peaks = pybedtools.BedTool(peaks_str, from_string=True).sort()
genes = pybedtools.BedTool(genes_str, from_string=True).sort()
with open('test_genome.txt', 'w') as fh:
    fh.write(genome_str)

print('=== Nearest gene per peak: -D b (sign by gene strand), -io, -t first ===')
near = peaks.closest(genes, D='b', io=True, t='first')
for iv in near:
    f = iv.fields
    dist = int(f[-1])
    if dist == -1:
        continue
    tag = 'DISTAL -> validate with ABC/PCHi-C/eQTL' if abs(dist) > DISTAL_FLAG else 'promoter-scale'
    print(f'  {f[3]} -> {f[9]} signed_dist={dist} ({tag})')

print('=== Candidate genes within a TAD-scale window (honest candidate set, counted) ===')
counts = peaks.window(genes, w=DISTAL_FLAG, c=True)
for iv in counts:
    print(f'  {iv.fields[3]}: {iv.fields[-1]} candidate gene(s) within {DISTAL_FLAG} bp')

print('=== Strand-aware promoters from TSS (collapse to TSS FIRST, then slop -s) ===')
tss = genes.each(lambda f: pybedtools.create_interval_from_list([f[0], str(f.start) if f.strand == '+' else str(f.end - 1), str(f.start + 1) if f.strand == '+' else str(f.end), f.name, f.score, f.strand])).saveas()
promoters = tss.slop(g='test_genome.txt', s=True, l=2000, r=200)
for iv in promoters:
    print(f'  {iv.name} ({iv.strand}): promoter {iv.chrom}:{iv.start}-{iv.end} width={iv.end - iv.start}')

import os
os.remove('test_genome.txt')
pybedtools.cleanup()
print('=== Done ===')
