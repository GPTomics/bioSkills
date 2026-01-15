#!/usr/bin/env python3
'''Count alleles at specific positions using pysam pileup'''

import pysam
import sys
from collections import Counter

def allele_counts(bam_path, chrom, pos, min_mapq=20, min_baseq=20):
    counts = Counter()

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for pileup_column in bam.pileup(chrom, pos, pos + 1,
                                         truncate=True,
                                         min_mapping_quality=min_mapq,
                                         min_base_quality=min_baseq):
            if pileup_column.pos != pos:
                continue

            for pileup_read in pileup_column.pileups:
                if pileup_read.is_del:
                    counts['DEL'] += 1
                elif pileup_read.is_refskip:
                    continue
                else:
                    qpos = pileup_read.query_position
                    base = pileup_read.alignment.query_sequence[qpos].upper()
                    counts[base] += 1

    return dict(counts)

def main():
    if len(sys.argv) < 3:
        print('Usage: allele_counts.py <input.bam> <region>')
        print('Example: allele_counts.py sample.bam chr1:1000000')
        sys.exit(1)

    bam_path = sys.argv[1]
    region = sys.argv[2]

    chrom, pos = region.split(':')
    pos = int(pos) - 1  # Convert to 0-based

    counts = allele_counts(bam_path, chrom, pos)
    total = sum(counts.values())

    print(f'Position: {chrom}:{pos+1}')
    print(f'Total depth: {total}')
    print('Allele counts:')

    for base, count in sorted(counts.items(), key=lambda x: -x[1]):
        freq = count / total * 100 if total > 0 else 0
        print(f'  {base}: {count} ({freq:.1f}%)')

if __name__ == '__main__':
    main()
