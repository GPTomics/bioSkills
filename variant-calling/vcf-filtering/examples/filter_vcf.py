#!/usr/bin/env python3
'''Filter VCF by quality criteria using cyvcf2'''

from cyvcf2 import VCF, Writer
import sys

def filter_vcf(input_path, output_path, min_qual=30, min_dp=10):
    vcf = VCF(input_path)
    writer = Writer(output_path, vcf)

    total = 0
    passed = 0

    for variant in vcf:
        total += 1
        qual = variant.QUAL or 0
        dp = variant.INFO.get('DP', 0)

        if qual >= min_qual and dp >= min_dp:
            writer.write_record(variant)
            passed += 1

    writer.close()
    vcf.close()

    print(f'Total variants: {total}')
    print(f'Passed filters: {passed}')
    print(f'Filtered out: {total - passed}')
    if total > 0:
        print(f'Pass rate: {passed/total*100:.1f}%')

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage: filter_vcf.py <input.vcf.gz> <output.vcf> [min_qual] [min_dp]')
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    min_qual = int(sys.argv[3]) if len(sys.argv) > 3 else 30
    min_dp = int(sys.argv[4]) if len(sys.argv) > 4 else 10

    filter_vcf(input_path, output_path, min_qual, min_dp)
