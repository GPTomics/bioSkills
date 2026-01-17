import pysam
import numpy as np
import argparse
from collections import Counter

def validate_bam(bam_file, sample_size=100000):
    bam = pysam.AlignmentFile(bam_file, 'rb')

    mapped, unmapped = 0, 0
    proper_pair, paired = 0, 0
    forward, reverse = 0, 0
    insert_sizes = []
    mapqs = []

    for i, read in enumerate(bam.fetch()):
        if i >= sample_size:
            break

        if read.is_unmapped:
            unmapped += 1
        else:
            mapped += 1
            forward += 1 if not read.is_reverse else 0
            reverse += 1 if read.is_reverse else 0
            mapqs.append(read.mapping_quality)

            if read.is_paired:
                paired += 1
                if read.is_proper_pair:
                    proper_pair += 1
                    if read.template_length > 0:
                        insert_sizes.append(read.template_length)

    bam.close()

    total = mapped + unmapped
    print(f'=== Alignment Validation (sampled {sample_size} reads) ===\n')

    print('--- Mapping ---')
    print(f'Mapped: {mapped} ({100*mapped/total:.1f}%)')
    print(f'Unmapped: {unmapped} ({100*unmapped/total:.1f}%)')

    print('\n--- Pairing ---')
    if paired > 0:
        print(f'Properly paired: {proper_pair} ({100*proper_pair/paired:.1f}%)')
    else:
        print('No paired reads')

    print('\n--- Insert Size ---')
    if insert_sizes:
        print(f'Median: {np.median(insert_sizes):.0f} bp')
        print(f'Mean: {np.mean(insert_sizes):.0f} bp')
        print(f'Std: {np.std(insert_sizes):.0f} bp')

    print('\n--- Strand Balance ---')
    print(f'Forward: {forward}, Reverse: {reverse}')
    print(f'Ratio: {forward/(forward+reverse):.3f}')

    print('\n--- MAPQ ---')
    print(f'Mean MAPQ: {np.mean(mapqs):.1f}')
    high_qual = sum(1 for m in mapqs if m >= 30)
    print(f'MAPQ >= 30: {100*high_qual/len(mapqs):.1f}%')

    print('\n--- Quality Summary ---')
    issues = []
    if mapped/total < 0.9:
        issues.append('Low mapping rate')
    if paired > 0 and proper_pair/paired < 0.8:
        issues.append('Low proper pairing')
    strand_ratio = forward/(forward+reverse)
    if strand_ratio < 0.45 or strand_ratio > 0.55:
        issues.append('Strand imbalance')
    if np.mean(mapqs) < 30:
        issues.append('Low mean MAPQ')

    if issues:
        print('WARNINGS:', ', '.join(issues))
    else:
        print('All metrics within normal range')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Validate BAM alignment quality')
    parser.add_argument('bam', help='Input BAM file')
    parser.add_argument('-n', '--sample-size', type=int, default=100000, help='Number of reads to sample')
    args = parser.parse_args()

    validate_bam(args.bam, args.sample_size)
