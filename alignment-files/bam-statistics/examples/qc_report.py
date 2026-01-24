#!/usr/bin/env python3
'''Generate QC report from BAM file'''

import pysam
import sys

def qc_report(bam_path):
    stats = {
        'total': 0, 'mapped': 0, 'paired': 0, 'proper_pair': 0,
        'duplicate': 0, 'secondary': 0, 'supplementary': 0
    }
    insert_sizes = []

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam:
            stats['total'] += 1
            if not read.is_unmapped:
                stats['mapped'] += 1
            if read.is_paired:
                stats['paired'] += 1
            if read.is_proper_pair:
                stats['proper_pair'] += 1
                if read.is_read1 and 0 < read.template_length < 1000:
                    insert_sizes.append(read.template_length)
            if read.is_duplicate:
                stats['duplicate'] += 1
            if read.is_secondary:
                stats['secondary'] += 1
            if read.is_supplementary:
                stats['supplementary'] += 1

    print(f'=== QC Report: {bam_path} ===\n')
    print(f'Total reads:       {stats["total"]:,}')
    print(f'Mapped:            {stats["mapped"]:,} ({stats["mapped"]/stats["total"]*100:.1f}%)')
    print(f'Properly paired:   {stats["proper_pair"]:,} ({stats["proper_pair"]/stats["paired"]*100:.1f}%)')
    print(f'Duplicates:        {stats["duplicate"]:,} ({stats["duplicate"]/stats["total"]*100:.1f}%)')
    print(f'Secondary:         {stats["secondary"]:,}')
    print(f'Supplementary:     {stats["supplementary"]:,}')

    if insert_sizes:
        mean_insert = sum(insert_sizes) / len(insert_sizes)
        insert_sizes.sort()
        median_insert = insert_sizes[len(insert_sizes) // 2]
        print(f'\nInsert size (mean): {mean_insert:.0f}')
        print(f'Insert size (median): {median_insert}')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: qc_report.py <input.bam>')
        sys.exit(1)

    qc_report(sys.argv[1])
