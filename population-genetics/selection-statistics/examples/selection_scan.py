#!/usr/bin/env python3
'''Genome-wide selection scan with Fst and Tajima's D.'''

import allel
import numpy as np
import matplotlib.pyplot as plt

vcf_file = 'data.vcf.gz'
pop1_samples = ['sample1', 'sample2', 'sample3', 'sample4', 'sample5']
pop2_samples = ['sample6', 'sample7', 'sample8', 'sample9', 'sample10']

print('=== Loading Data ===')
callset = allel.read_vcf(vcf_file)
gt = allel.GenotypeArray(callset['calldata/GT'])
pos = callset['variants/POS']
samples = list(callset['samples'])

pop1_idx = [samples.index(s) for s in pop1_samples if s in samples]
pop2_idx = [samples.index(s) for s in pop2_samples if s in samples]

print(f'Pop1 samples: {len(pop1_idx)}')
print(f'Pop2 samples: {len(pop2_idx)}')

print('\n=== Filtering ===')
ac = gt.count_alleles()
flt = ac.is_segregating() & (ac.max_allele() == 1)
gt = gt.compress(flt, axis=0)
pos = pos[flt]
print(f'Variants after filtering: {gt.n_variants}')

print('\n=== Calculating Fst ===')
subpops = {'pop1': pop1_idx, 'pop2': pop2_idx}
ac_subpops = gt.count_alleles_subpops(subpops)

fst, fst_windows, _ = allel.windowed_hudson_fst(
    pos, ac_subpops['pop1'], ac_subpops['pop2'],
    size=100000, step=50000)

print(f'Mean genome-wide Fst: {np.nanmean(fst):.4f}')

print("\n=== Calculating Tajima's D ===")
ac_all = gt.count_alleles()
tajD, tajD_windows, _ = allel.windowed_tajima_d(pos, ac_all, size=100000, step=50000)
print(f"Mean Tajima's D: {np.nanmean(tajD):.4f}")

print('\n=== Plotting ===')
fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

axes[0].plot(fst_windows[:, 0] / 1e6, fst, 'b-', linewidth=0.5)
axes[0].axhline(np.nanmean(fst), color='r', linestyle='--', alpha=0.5)
axes[0].set_ylabel('Fst')
axes[0].set_title('Population Differentiation')

axes[1].plot(tajD_windows[:, 0] / 1e6, tajD, 'g-', linewidth=0.5)
axes[1].axhline(0, color='r', linestyle='--', alpha=0.5)
axes[1].set_ylabel("Tajima's D")
axes[1].set_xlabel('Position (Mb)')
axes[1].set_title('Departures from Neutrality')

plt.tight_layout()
plt.savefig('selection_scan.png', dpi=150)
print('Plot saved: selection_scan.png')

print('\n=== Top Fst Windows ===')
top_fst_idx = np.argsort(fst)[-10:][::-1]
for i in top_fst_idx:
    if not np.isnan(fst[i]):
        print(f'{fst_windows[i, 0]:,}-{fst_windows[i, 1]:,}: Fst={fst[i]:.4f}')

print("\n=== Extreme Tajima's D Windows ===")
low_tajD_idx = np.argsort(tajD)[:5]
for i in low_tajD_idx:
    if not np.isnan(tajD[i]):
        print(f'{tajD_windows[i, 0]:,}-{tajD_windows[i, 1]:,}: D={tajD[i]:.4f}')

print('\n=== Done ===')
