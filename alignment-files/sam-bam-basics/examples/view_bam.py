#!/usr/bin/env python3
'''View BAM file contents using pysam'''

import pysam
import sys

def view_bam(bam_path, limit=10):
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        print(f'References: {bam.nreferences}')
        print(f'Mapped: {bam.mapped}')
        print(f'Unmapped: {bam.unmapped}')
        print()

        for i, read in enumerate(bam):
            if i >= limit:
                break
            strand = '-' if read.is_reverse else '+'
            print(f'{read.query_name}\t{read.reference_name}:{read.reference_start}\t{strand}\t{read.cigarstring}')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: view_bam.py <input.bam> [limit]')
        sys.exit(1)

    bam_path = sys.argv[1]
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    view_bam(bam_path, limit)
