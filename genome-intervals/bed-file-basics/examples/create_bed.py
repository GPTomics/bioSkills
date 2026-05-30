#!/usr/bin/env python3
'''Create BED files and round-trip a 1-based -> 0-based coordinate conversion.'''
# Reference: pybedtools 0.10+, bedtools 2.31+, pandas 2.2+ | Verify API if version differs

import pybedtools
import pandas as pd

intervals = [('chr1', 1000, 2000, 'gene1', 500, '+'), ('chr1', 3000, 4000, 'gene2', 600, '-'), ('chr2', 5000, 6000, 'gene3', 450, '+')]
bed = pybedtools.BedTool(intervals).saveas('from_list.bed')
print(f'from_list.bed: {bed.count()} intervals')

df = pd.DataFrame({'chrom': ['chr1', 'chr1', 'chr2'], 'start': [100, 500, 200], 'end': [200, 700, 400], 'name': ['r1', 'r2', 'r3'], 'score': [100, 200, 150], 'strand': ['+', '-', '+']})
bed_df = pybedtools.BedTool.from_dataframe(df).saveas('from_dataframe.bed')
print(f'from_dataframe.bed: {bed_df.count()} intervals')

# Convert a 1-based GTF/VCF position to BED: subtract 1 from start, leave end unchanged.
# A 1 bp feature at 1-based position 6 (GTF 'chr1 6 6') is BED 'chr1 5 6' -- length stays 1.
gtf_start_1based, gtf_end_1based = 6, 6
bed_start, bed_end = gtf_start_1based - 1, gtf_end_1based
assert bed_end - bed_start == 1, 'round-trip failed: a 1 bp feature must stay 1 bp'
print(f'1-based [{gtf_start_1based},{gtf_end_1based}] -> BED [{bed_start},{bed_end}); length {bed_end - bed_start}')

# len(interval) reports end - start (no +1); BED has no fencepost.
for iv in bed_df:
    print(iv.chrom, iv.start, iv.end, len(iv))

pybedtools.cleanup()
