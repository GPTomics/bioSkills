#!/usr/bin/env python3
'''Basic population genetics analysis with scikit-allel.'''

import allel
import numpy as np
import matplotlib.pyplot as plt

vcf_file = 'data.vcf.gz'

print('=== Loading VCF ===')
callset = allel.read_vcf(vcf_file,
    fields=['samples', 'calldata/GT', 'variants/POS', 'variants/CHROM'])

gt = allel.GenotypeArray(callset['calldata/GT'])
pos = callset['variants/POS']
samples = callset['samples']

print(f'Samples: {len(samples)}')
print(f'Variants: {gt.n_variants}')

print('\n=== Quality Filtering ===')
ac = gt.count_alleles()
is_seg = ac.is_segregating()
is_biallelic = ac.max_allele() == 1

flt = is_seg & is_biallelic
gt = gt.compress(flt, axis=0)
pos = pos[flt]
ac = gt.count_alleles()

print(f'Variants after filtering: {gt.n_variants}')

print('\n=== Diversity Statistics ===')
pi = allel.sequence_diversity(pos, ac)
theta_w = allel.watterson_theta(pos, ac)
ho = allel.heterozygosity_observed(gt).mean()

print(f'Nucleotide diversity (Pi): {pi:.6f}')
print(f"Watterson's theta: {theta_w:.6f}")
print(f'Mean observed heterozygosity: {ho:.4f}')

print('\n=== PCA ===')
gn = gt.to_n_alt(fill=-1)
gn = np.where(gn < 0, 0, gn)  # Impute missing

coords, model = allel.pca(gn, n_components=10, scaler='patterson')

plt.figure(figsize=(8, 6))
plt.scatter(coords[:, 0], coords[:, 1], s=20)
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('PCA of Genotype Data')
plt.savefig('pca_result.png', dpi=150)
print('PCA plot saved: pca_result.png')

print('\n=== Site Frequency Spectrum ===')
sfs = allel.sfs(ac[:, 1])
sfs_folded = allel.sfs_folded(ac)

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].bar(range(len(sfs)), sfs)
axes[0].set_xlabel('Derived allele count')
axes[0].set_ylabel('Frequency')
axes[0].set_title('Unfolded SFS')

axes[1].bar(range(len(sfs_folded)), sfs_folded)
axes[1].set_xlabel('Minor allele count')
axes[1].set_ylabel('Frequency')
axes[1].set_title('Folded SFS')

plt.tight_layout()
plt.savefig('sfs_result.png', dpi=150)
print('SFS plot saved: sfs_result.png')

print('\n=== Done ===')
