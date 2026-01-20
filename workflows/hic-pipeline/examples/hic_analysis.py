# Hi-C analysis: compartments, TADs, and loops

import cooler
import cooltools
import cooltools.lib.plotting
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Configuration
mcool_path = 'sample.mcool'
output_dir = 'hic_analysis'
os.makedirs(output_dir, exist_ok=True)
os.makedirs(f'{output_dir}/plots', exist_ok=True)

# === Compartments (100kb resolution) ===
print('=== Compartment Analysis ===')
clr_100k = cooler.Cooler(f'{mcool_path}::resolutions/100000')
print(f'Matrix shape: {clr_100k.shape}')

# Calculate expected for O/E
expected_100k = cooltools.expected_cis(clr_100k, nproc=4)

# Eigenvector decomposition
eig_values, eig_vectors = cooltools.eigs_cis(
    clr_100k,
    n_eigs=3,
    nproc=4
)

# Assign A/B compartments
eig_vectors['compartment'] = np.where(eig_vectors['E1'] > 0, 'A', 'B')
eig_vectors.to_csv(f'{output_dir}/compartments.tsv', sep='\t', index=False)

print(f'A compartment bins: {(eig_vectors["compartment"] == "A").sum()}')
print(f'B compartment bins: {(eig_vectors["compartment"] == "B").sum()}')

# Plot compartment track for chr1
chr1_eig = eig_vectors[eig_vectors['chrom'] == 'chr1']
fig, ax = plt.subplots(figsize=(12, 3))
ax.fill_between(range(len(chr1_eig)), chr1_eig['E1'],
                where=chr1_eig['E1'] > 0, color='red', alpha=0.7, label='A')
ax.fill_between(range(len(chr1_eig)), chr1_eig['E1'],
                where=chr1_eig['E1'] < 0, color='blue', alpha=0.7, label='B')
ax.axhline(0, color='black', linewidth=0.5)
ax.set_xlabel('Position (100kb bins)')
ax.set_ylabel('E1 (compartment score)')
ax.set_title('A/B Compartments - chr1')
ax.legend()
plt.tight_layout()
plt.savefig(f'{output_dir}/plots/compartments_chr1.pdf')
plt.close()

# === TAD Detection (10kb resolution) ===
print('\n=== TAD Detection ===')
clr_10k = cooler.Cooler(f'{mcool_path}::resolutions/10000')

# Calculate insulation score
insulation = cooltools.insulation(
    clr_10k,
    window_bp=[100000, 200000, 500000],
    nproc=4
)
insulation.to_csv(f'{output_dir}/insulation_scores.tsv', sep='\t', index=False)

# Find boundaries (minima in insulation score)
boundaries = insulation.copy()
boundaries['is_boundary'] = False

# Simple boundary detection
for chrom in boundaries['chrom'].unique():
    if pd.isna(chrom):
        continue
    chrom_data = boundaries[boundaries['chrom'] == chrom]
    insul_col = 'log2_insulation_score_100000'
    if insul_col in chrom_data.columns:
        # Find local minima
        values = chrom_data[insul_col].values
        for i in range(1, len(values) - 1):
            if values[i] < values[i-1] and values[i] < values[i+1]:
                if values[i] < -0.1:  # Threshold
                    boundaries.loc[chrom_data.index[i], 'is_boundary'] = True

tad_boundaries = boundaries[boundaries['is_boundary']]
tad_boundaries.to_csv(f'{output_dir}/tad_boundaries.tsv', sep='\t', index=False)
print(f'TAD boundaries found: {len(tad_boundaries)}')

# Plot insulation for chr1
chr1_ins = insulation[insulation['chrom'] == 'chr1'].iloc[:500]
fig, ax = plt.subplots(figsize=(12, 3))
ax.plot(chr1_ins['log2_insulation_score_100000'], color='black')
ax.axhline(0, color='gray', linestyle='--', linewidth=0.5)
ax.set_xlabel('Position (10kb bins)')
ax.set_ylabel('log2(insulation score)')
ax.set_title('Insulation Score - chr1')
plt.tight_layout()
plt.savefig(f'{output_dir}/plots/insulation_chr1.pdf')
plt.close()

# === Loop Calling (10kb resolution) ===
print('\n=== Loop Calling ===')
expected_10k = cooltools.expected_cis(clr_10k, nproc=4)

# Call loops
loops = cooltools.dots(
    clr_10k,
    expected_10k,
    max_loci_separation=2000000,
    nproc=4
)
loops.to_csv(f'{output_dir}/loops.tsv', sep='\t', index=False)
print(f'Loops called: {len(loops)}')

# === Summary ===
print('\n=== Analysis Complete ===')
print(f'Results saved to: {output_dir}/')
print(f'  - compartments.tsv: A/B compartment assignments')
print(f'  - insulation_scores.tsv: TAD insulation')
print(f'  - tad_boundaries.tsv: TAD boundary positions')
print(f'  - loops.tsv: Chromatin loops')
print(f'  - plots/: Visualization PDFs')
