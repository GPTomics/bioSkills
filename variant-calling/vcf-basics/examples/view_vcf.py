#!/usr/bin/env python3
'''View VCF file contents using cyvcf2'''

from cyvcf2 import VCF
import sys

def view_vcf(vcf_path, limit=10):
    vcf = VCF(vcf_path)

    print(f'Samples: {", ".join(vcf.samples)}')
    print(f'Contigs: {len(vcf.seqnames)}')
    print()

    for i, variant in enumerate(vcf):
        if i >= limit:
            break

        alt = ','.join(variant.ALT) if variant.ALT else '.'
        qual = f'{variant.QUAL:.1f}' if variant.QUAL else '.'
        filt = variant.FILTER if variant.FILTER else 'PASS'

        print(f'{variant.CHROM}:{variant.POS}\t{variant.REF}>{alt}\t'
              f'QUAL={qual}\tFILTER={filt}\tTYPE={variant.var_type}')

    vcf.close()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: view_vcf.py <input.vcf.gz> [limit]')
        sys.exit(1)

    vcf_path = sys.argv[1]
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    view_vcf(vcf_path, limit)
