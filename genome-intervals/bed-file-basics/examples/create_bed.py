#!/usr/bin/env python3
'''Create BED files from various sources using pybedtools.'''

import pybedtools
import pandas as pd

# Method 1: From string
bed_string = '''chr1\t100\t200\tpeak1\t100\t+
chr1\t300\t400\tpeak2\t200\t-
chr2\t500\t600\tpeak3\t150\t+'''
bed1 = pybedtools.BedTool(bed_string, from_string=True)
bed1.saveas('from_string.bed')
print(f'Created from_string.bed with {bed1.count()} intervals')

# Method 2: From list of tuples
intervals = [
    ('chr1', 1000, 2000, 'gene1', 500, '+'),
    ('chr1', 3000, 4000, 'gene2', 600, '-'),
    ('chr2', 5000, 6000, 'gene3', 450, '+'),
]
bed2 = pybedtools.BedTool(intervals)
bed2.saveas('from_list.bed')
print(f'Created from_list.bed with {bed2.count()} intervals')

# Method 3: From pandas DataFrame
df = pd.DataFrame({
    'chrom': ['chr1', 'chr1', 'chr2', 'chr3'],
    'start': [100, 500, 200, 1000],
    'end': [200, 700, 400, 2000],
    'name': ['region1', 'region2', 'region3', 'region4'],
    'score': [100, 200, 150, 300],
    'strand': ['+', '-', '+', '-']
})
bed3 = pybedtools.BedTool.from_dataframe(df)
bed3.saveas('from_dataframe.bed')
print(f'Created from_dataframe.bed with {bed3.count()} intervals')

# Method 4: Create BED3 minimal format
bed3_string = '''chr1\t0\t1000
chr1\t2000\t3000
chr2\t0\t5000'''
bed4 = pybedtools.BedTool(bed3_string, from_string=True)
bed4.saveas('minimal.bed')
print(f'Created minimal.bed with {bed4.count()} intervals')

# Read back and convert to DataFrame
print('\nContents of from_dataframe.bed:')
bed_loaded = pybedtools.BedTool('from_dataframe.bed')
df_loaded = bed_loaded.to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
print(df_loaded)

# Cleanup temp files
pybedtools.cleanup()
