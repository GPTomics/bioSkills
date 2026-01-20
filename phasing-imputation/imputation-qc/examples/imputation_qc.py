#!/usr/bin/env python3
'''QC metrics for imputation results'''

from cyvcf2 import VCF
import numpy as np
import matplotlib.pyplot as plt

def imputation_qc(vcf_path, output_prefix='imputation_qc'):
    '''Calculate imputation QC metrics'''
    vcf = VCF(vcf_path)

    r2_values = []
    maf_values = []
    n_variants = 0

    for v in vcf:
        n_variants += 1
        r2 = v.INFO.get('DR2', v.INFO.get('R2', None))
        af = v.INFO.get('AF', None)

        if r2 is not None:
            r2_values.append(r2)
        if af is not None:
            maf = min(af, 1 - af) if isinstance(af, float) else min(af[0], 1 - af[0])
            maf_values.append(maf)

    vcf.close()

    print(f'Total variants: {n_variants}')
    print(f'Mean R2: {np.mean(r2_values):.3f}')
    print(f'Variants with R2 > 0.8: {sum(r > 0.8 for r in r2_values)} ({100*sum(r > 0.8 for r in r2_values)/len(r2_values):.1f}%)')
    print(f'Variants with R2 > 0.3: {sum(r > 0.3 for r in r2_values)} ({100*sum(r > 0.3 for r in r2_values)/len(r2_values):.1f}%)')

    # Plot R2 distribution
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].hist(r2_values, bins=50, edgecolor='black')
    axes[0].axvline(0.8, color='red', linestyle='--', label='R2=0.8')
    axes[0].set_xlabel('Imputation R2')
    axes[0].set_ylabel('Count')
    axes[0].legend()

    axes[1].scatter(maf_values[:10000], r2_values[:10000], alpha=0.3, s=5)
    axes[1].set_xlabel('MAF')
    axes[1].set_ylabel('R2')
    axes[1].set_title('R2 vs MAF')

    plt.tight_layout()
    plt.savefig(f'{output_prefix}.png', dpi=150)
    print(f'Saved: {output_prefix}.png')

if __name__ == '__main__':
    import sys
    vcf_path = sys.argv[1] if len(sys.argv) > 1 else 'imputed.vcf.gz'
    imputation_qc(vcf_path)
