#!/usr/bin/env python3
'''Generate VCF statistics using cyvcf2'''

from cyvcf2 import VCF
import sys

def calculate_stats(vcf_path):
    vcf = VCF(vcf_path)

    stats = {
        'total': 0, 'snps': 0, 'indels': 0, 'other': 0,
        'pass': 0, 'filtered': 0,
        'transitions': 0, 'transversions': 0,
        'qual_sum': 0, 'qual_count': 0
    }

    ti_pairs = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}

    for variant in vcf:
        stats['total'] += 1

        if variant.is_snp:
            stats['snps'] += 1
            ref, alt = variant.REF, variant.ALT[0]
            if (ref, alt) in ti_pairs:
                stats['transitions'] += 1
            else:
                stats['transversions'] += 1
        elif variant.is_indel:
            stats['indels'] += 1
        else:
            stats['other'] += 1

        if variant.FILTER is None:
            stats['pass'] += 1
        else:
            stats['filtered'] += 1

        if variant.QUAL:
            stats['qual_sum'] += variant.QUAL
            stats['qual_count'] += 1

    vcf.close()
    return stats

def print_stats(stats):
    print('=== VCF Statistics ===')
    print(f'Total variants:    {stats["total"]}')
    print(f'  SNPs:            {stats["snps"]}')
    print(f'  Indels:          {stats["indels"]}')
    print(f'  Other:           {stats["other"]}')
    print()
    print(f'PASS variants:     {stats["pass"]}')
    print(f'Filtered:          {stats["filtered"]}')
    print()
    if stats['transversions'] > 0:
        tstv = stats['transitions'] / stats['transversions']
        print(f'Transitions:       {stats["transitions"]}')
        print(f'Transversions:     {stats["transversions"]}')
        print(f'Ti/Tv ratio:       {tstv:.2f}')
    print()
    if stats['qual_count'] > 0:
        mean_qual = stats['qual_sum'] / stats['qual_count']
        print(f'Mean QUAL:         {mean_qual:.1f}')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: vcf_stats.py <input.vcf.gz>')
        sys.exit(1)

    vcf_path = sys.argv[1]
    stats = calculate_stats(vcf_path)
    print_stats(stats)
