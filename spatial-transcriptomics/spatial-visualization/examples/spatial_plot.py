'''Fork-aware, honest spatial visualization on a built-in Visium section.

Demonstrates spot/capture plotting with the dataset scalefactors, a perceptually
uniform colormap, disclosed clipping, and NO interpolation of the sparse field.
Figures are written to a temp dir so the repo stays clean.
'''
# Reference: squidpy 1.4+, scanpy 1.10+, anndata 0.10+, matplotlib 3.8+ | Verify API if version differs

import tempfile
import os
import scanpy as sc
import squidpy as sq

adata = sq.datasets.visium_hne_adata()                 # Visium H&E section: image + scalefactors in uns['spatial']
print(f'Loaded: {adata.n_obs} spots, {adata.n_vars} genes')

sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None)
gene = adata.var_names[0]                               # any present gene; imaging panels would constrain this set

outdir = tempfile.mkdtemp(prefix='spatial_viz_')

# Spot/capture fork: sc.pl.spatial reads scalefactors so spots align to the histology image.
# A spot is a 1-10-cell MIXTURE, not a cell -- cluster labels here mark regions, not cell types.
sc.pl.spatial(adata, color='cluster', img_key='hires', alpha_img=0.6, show=False,
              save=None, title='Clusters (regions, not cell types)')

# Honest continuous field: color measured spots only (no KDE/contour), perceptually uniform map,
# disclosed p99 clip. Oversized markers would fake coverage, so keep the default spot scaling.
sc.pl.spatial(adata, color=['total_counts', gene], cmap='viridis', vmin=0, vmax='p99',
              ncols=2, show=False)

# Coordinate-frame transform: hires-image pixel coords = spatial coords * tissue_hires_scalef.
library_id = list(adata.uns['spatial'].keys())[0]
scalef = adata.uns['spatial'][library_id]['scalefactors']['tissue_hires_scalef']
coords_px = adata.obsm['spatial'] * scalef
print(f'First spot microns->pixels: {adata.obsm["spatial"][0]} -> {coords_px[0]}')

fig_path = os.path.join(outdir, 'spatial_overview.png')
import matplotlib.pyplot as plt
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
print(f'Wrote figure to temp dir: {fig_path}')
