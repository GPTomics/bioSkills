'''Call A/B compartments from Hi-C data'''

import cooler
import cooltools
import bioframe
import numpy as np
import matplotlib.pyplot as plt

clr = cooler.Cooler('matrix.mcool::resolutions/100000')
print(f'Loaded at {clr.binsize}bp resolution')

view_df = bioframe.make_viewframe(clr.chromsizes)

print('Computing expected values...')
expected = cooltools.expected_cis(clr, view_df=view_df, ignore_diags=2)

print('Computing eigenvectors...')
eigenvalues, eigenvectors = cooltools.eigs_cis(clr, view_df=view_df, n_eigs=1)

eigenvectors['compartment'] = np.where(eigenvectors['E1'] > 0, 'A', 'B')
print(f'\nCompartment counts:')
print(eigenvectors['compartment'].value_counts())

eigenvectors[['chrom', 'start', 'end', 'E1', 'compartment']].to_csv(
    'compartments.bed', sep='\t', index=False, header=False
)
print('\nSaved to compartments.bed')

print('\nComputing saddle plot...')
saddle_data = cooltools.saddle(
    clr, expected=expected, eigenvector_track=eigenvectors,
    view_df=view_df, n_bins=50, vrange=(-0.5, 0.5)
)

saddle_agg = np.nanmean(saddle_data[0], axis=0)
plt.figure(figsize=(6, 6))
plt.imshow(saddle_agg, cmap='coolwarm', vmin=-1, vmax=1)
plt.colorbar(label='log2(O/E)')
plt.title('Saddle plot')
plt.savefig('saddle_plot.png', dpi=150)
print('Saved saddle_plot.png')
