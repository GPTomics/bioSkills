#!/usr/bin/env python3
'''Filter BAM file by various criteria'''

import pysam
import argparse

def filter_bam(input_bam, output_bam, min_mapq=0, remove_duplicates=False,
               primary_only=False, proper_pair=False, region=None):

    with pysam.AlignmentFile(input_bam, 'rb') as infile:
        with pysam.AlignmentFile(output_bam, 'wb', header=infile.header) as outfile:
            if region:
                chrom, coords = region.split(':')
                start, end = map(int, coords.split('-'))
                iterator = infile.fetch(chrom, start, end)
            else:
                iterator = infile

            kept = 0
            removed = 0

            for read in iterator:
                if read.is_unmapped:
                    removed += 1
                    continue
                if read.mapping_quality < min_mapq:
                    removed += 1
                    continue
                if remove_duplicates and read.is_duplicate:
                    removed += 1
                    continue
                if primary_only and (read.is_secondary or read.is_supplementary):
                    removed += 1
                    continue
                if proper_pair and not read.is_proper_pair:
                    removed += 1
                    continue

                outfile.write(read)
                kept += 1

    print(f'Kept: {kept:,}')
    print(f'Removed: {removed:,}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter BAM file')
    parser.add_argument('input', help='Input BAM file')
    parser.add_argument('output', help='Output BAM file')
    parser.add_argument('-q', '--min-mapq', type=int, default=0, help='Minimum MAPQ')
    parser.add_argument('-d', '--remove-duplicates', action='store_true', help='Remove duplicates')
    parser.add_argument('-p', '--primary-only', action='store_true', help='Primary alignments only')
    parser.add_argument('-P', '--proper-pair', action='store_true', help='Properly paired only')
    parser.add_argument('-r', '--region', help='Region (chr:start-end)')

    args = parser.parse_args()

    filter_bam(args.input, args.output, args.min_mapq, args.remove_duplicates,
               args.primary_only, args.proper_pair, args.region)

    print('\nIndexing output...')
    pysam.index(args.output)
