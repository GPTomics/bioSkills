'''Call GC-phased A/B compartments per chromosome arm and quantify compartment strength via a saddle plot.'''
# Reference: cooler 0.10+, cooltools 0.7+, bioframe 0.7+, matplotlib 3.8+ | Verify API if version differs

import cooler
import cooltools
import bioframe
import numpy as np
import matplotlib.pyplot as plt

COMPARTMENT_RES = 100_000   # 100kb: compartments are chromosome-scale; finer bins mix in TADs/loops and sparsity noise
N_EIGS = 3                  # request >=3 so the arm-gradient-vs-compartment problem is visible in the eigenvalues
N_GROUPS = 38               # quantile groups for the saddle; ~30-50 is conventional (cooltools tutorial)
Q_LO, Q_HI = 0.025, 0.975   # trim the extreme 2.5% tails before digitizing to resist outlier bins
GENOME = 'hg38'

clr = cooler.Cooler(f'matrix.mcool::/resolutions/{COMPARTMENT_RES}')
print(f'Loaded at {clr.binsize}bp resolution')

cens = bioframe.fetch_centromeres(GENOME)
arms = bioframe.make_chromarms(clr.chromsizes, cens)            # per-arm view removes the centromere gradient
arms = arms[arms.chrom.isin(clr.chromnames)].reset_index(drop=True)

genome = bioframe.load_fasta(f'{GENOME}.fa')                    # indexed FASTA (.fai must exist)
bins = clr.bins()[:][['chrom', 'start', 'end']]
gc = bioframe.frac_gc(bins, genome)                            # GC phasing track at the cooler's exact binning

print('Computing GC-phased eigenvectors per arm...')
eigvals, eigvecs = cooltools.eigs_cis(clr, gc, view_df=arms, n_eigs=N_EIGS, sort_metric='pearsonr')
eigvecs['compartment'] = ['A' if e > 0 else 'B' for e in eigvecs['E1']]
print(eigvecs['compartment'].value_counts())

corr = np.corrcoef(eigvecs['E1'].fillna(0), gc['GC'].fillna(0))[0, 1]
print(f'E1 vs GC correlation: {corr:.3f} (should be strongly positive; if weak, the compartment signal may be in E2/E3)')

eigvecs[['chrom', 'start', 'end', 'E1', 'compartment']].to_csv('compartments.bed', sep='\t', index=False, header=False)
print('Saved compartments.bed')

print('Computing saddle plot and strength...')
expected = cooltools.expected_cis(clr, view_df=arms)           # provides the 'balanced.avg' column saddle needs
track = eigvecs[['chrom', 'start', 'end', 'E1']]               # the SAME phased E1 used for the A/B call
interaction_sum, interaction_count = cooltools.saddle(clr, expected, track, 'cis', n_bins=N_GROUPS, qrange=(Q_LO, Q_HI), view_df=arms)

strength = cooltools.api.saddle.saddle_strength(interaction_sum, interaction_count)   # 1D array over corner extent; not exposed at top level
extent = N_GROUPS // 5   # read strength at the top/bottom ~20% of bins; keep this extent identical across compared samples
print(f'Compartment strength at extent {extent}: {strength[extent]:.3f}')

saddle = interaction_sum / interaction_count
plt.figure(figsize=(6, 6))
plt.imshow(np.log2(saddle), cmap='coolwarm', vmin=-1, vmax=1)
plt.colorbar(label='log2(O/E)')
plt.xlabel('E1 quantile (B -> A)')
plt.ylabel('E1 quantile (B -> A)')
plt.title('Saddle plot')
plt.savefig('saddle_plot.png', dpi=150)
print('Saved saddle_plot.png')
