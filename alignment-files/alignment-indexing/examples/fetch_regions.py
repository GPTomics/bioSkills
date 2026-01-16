#!/usr/bin/env python3
'''Fetch reads from specific regions using index'''

import pysam
import sys
from pathlib import Path

def ensure_indexed(bam_path):
    bam_path = Path(bam_path)
    if not (bam_path.with_suffix('.bam.bai').exists() or
            Path(str(bam_path) + '.bai').exists()):
        print(f'Indexing {bam_path}...')
        pysam.index(str(bam_path))

def fetch_region(bam_path, region):
    ensure_indexed(bam_path)

    chrom, coords = region.split(':')
    start, end = map(int, coords.split('-'))

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        count = 0
        for read in bam.fetch(chrom, start - 1, end):  # Convert to 0-based
            count += 1
            strand = '-' if read.is_reverse else '+'
            print(f'{read.query_name}\t{read.reference_start + 1}\t{strand}')

        print(f'\nTotal reads in {region}: {count}')

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage: fetch_regions.py <input.bam> <region>')
        print('Example: fetch_regions.py sample.bam chr1:1000000-2000000')
        sys.exit(1)

    fetch_region(sys.argv[1], sys.argv[2])
