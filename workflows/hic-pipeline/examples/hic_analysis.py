# Reference: cooler 0.10+, cooltools 0.7+, bioframe 0.7+, matplotlib 3.8+ | Verify API if version differs
# Hi-C analysis: compartments, TADs, and loops from a balanced multi-resolution cooler.

import cooler
import cooltools
import bioframe
import numpy as np
import matplotlib.pyplot as plt
import os

mcool_path = 'sample.mcool'
genome_fasta = 'reference.fa'
output_dir = 'hic_analysis'
os.makedirs(f'{output_dir}/plots', exist_ok=True)

COMPARTMENT_RES = 100_000   # 100kb: compartments are chromosome-scale; finer bins mix in TADs
TAD_RES = 10_000            # 10kb: domain/boundary scale
INSULATION_WINDOWS = [100_000, 200_000, 500_000]   # diamond windows; the window IS the scale dial
BOUNDARY_WINDOW = 200_000   # report boundaries called at the mid window
MAX_LOOP_SEP = 2_000_000    # dots beyond 2Mb are dominated by noise at typical depth

# === Compartments (100kb) -- sign-phased by GC, else A/B is arbitrary ===
print('=== Compartment Analysis ===')
clr_100k = cooler.Cooler(f'{mcool_path}::/resolutions/{COMPARTMENT_RES}')
gc = bioframe.frac_gc(clr_100k.bins()[:][['chrom', 'start', 'end']], bioframe.load_fasta(genome_fasta))
eig_values, eig_vectors = cooltools.eigs_cis(clr_100k, gc, n_eigs=3)
# Leave masked/unmappable bins (E1 = NaN from mad_max-dropped regions) unassigned; NaN > 0 is
# False, so a bare np.where would silently label every unmappable bin 'B'.
eig_vectors['compartment'] = np.where(eig_vectors['E1'].isna(), None, np.where(eig_vectors['E1'] > 0, 'A', 'B'))
eig_vectors.to_csv(f'{output_dir}/compartments.tsv', sep='\t', index=False)
print(f'A bins: {(eig_vectors["compartment"] == "A").sum()}  B bins: {(eig_vectors["compartment"] == "B").sum()}')

chr1_eig = eig_vectors[eig_vectors['chrom'] == 'chr1']
fig, ax = plt.subplots(figsize=(12, 3))
ax.fill_between(range(len(chr1_eig)), chr1_eig['E1'], where=chr1_eig['E1'] > 0, color='red', alpha=0.7, label='A')
ax.fill_between(range(len(chr1_eig)), chr1_eig['E1'], where=chr1_eig['E1'] < 0, color='blue', alpha=0.7, label='B')
ax.axhline(0, color='black', linewidth=0.5)
ax.set(xlabel='Position (100kb bins)', ylabel='E1', title='A/B Compartments - chr1')
ax.legend()
plt.tight_layout(); plt.savefig(f'{output_dir}/plots/compartments_chr1.pdf'); plt.close()

# === TAD boundaries (10kb) -- insulation() returns boundary columns directly ===
print('=== TAD Detection ===')
clr_10k = cooler.Cooler(f'{mcool_path}::/resolutions/{TAD_RES}')
ins = cooltools.insulation(clr_10k, window_bp=INSULATION_WINDOWS)
ins.to_csv(f'{output_dir}/insulation_scores.tsv', sep='\t', index=False)

is_col, strength_col = f'is_boundary_{BOUNDARY_WINDOW}', f'boundary_strength_{BOUNDARY_WINDOW}'
# is_boundary is NaN for bad/low-mappability bins; fillna(False) before masking or pandas raises on the NaN.
boundaries = ins[ins[is_col].fillna(False).astype(bool)][['chrom', 'start', 'end', strength_col]]
boundaries.to_csv(f'{output_dir}/tad_boundaries.tsv', sep='\t', index=False)
print(f'TAD boundaries ({BOUNDARY_WINDOW}bp window): {len(boundaries)}')

chr1_ins = ins[ins['chrom'] == 'chr1'].iloc[:500]
fig, ax = plt.subplots(figsize=(12, 3))
ax.plot(chr1_ins[f'log2_insulation_score_{BOUNDARY_WINDOW}'].values, color='black')
ax.axhline(0, color='gray', linestyle='--', linewidth=0.5)
ax.set(xlabel='Position (10kb bins)', ylabel='log2(insulation)', title='Insulation Score - chr1')
plt.tight_layout(); plt.savefig(f'{output_dir}/plots/insulation_chr1.pdf'); plt.close()

# === Loop calling (10kb) -- only honest on a deep map; else run APA on known anchors ===
print('=== Loop Calling ===')
expected_10k = cooltools.expected_cis(clr_10k, nproc=4)
loops = cooltools.dots(clr_10k, expected_10k, max_loci_separation=MAX_LOOP_SEP, nproc=4)
loops.to_csv(f'{output_dir}/loops.tsv', sep='\t', index=False)
print(f'Loops called: {len(loops)}')

print(f'Results saved to: {output_dir}/')
