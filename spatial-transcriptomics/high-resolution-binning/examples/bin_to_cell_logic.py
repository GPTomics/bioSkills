'''Demonstrate the bin-UP-into-cells logic on a synthetic 2um grid (no bin2cell/GPU needed).

Shows the two facts the skill turns on:
  1. assigning sub-cellular bins to nuclei and summing them reconstructs a cell (the Bin2cell idea);
  2. a coarse fixed grid re-creates the multi-cell mixture that reconstruction avoids.
The real workflow uses b2c.read_visium -> b2c.stardist -> b2c.insert_labels -> b2c.bin_to_cell;
this is a dependency-free illustration of the same aggregation, not a substitute for it.
'''
# Reference: bin2cell 0.3+ | Verify API if version differs
import numpy as np
import pandas as pd

rng = np.random.default_rng(0)

# two real cells, ~10um apart, each a nucleus center on a micron grid
nuclei = np.array([[20.0, 20.0], [34.0, 20.0]])   # cell 0 and cell 1, ~14um apart (touching at ~10um radii)
cell_radius_um = 7.0                              # nucleus-expansion territory; bins beyond all territories are intercellular

# distinct expression identity per cell over a tiny 4-gene panel
genes = ['Epcam', 'Cd3d', 'Col1a1', 'Pecam1']
cell_profiles = np.array([[8.0, 0.5, 0.2, 0.2],   # cell 0: epithelial
                          [0.3, 9.0, 0.2, 0.4]])  # cell 1: T cell

# lay down a 2um bin lattice over the field, each bin a sparse fragment of whichever cell owns it
step = 2.0
xs, ys = np.meshgrid(np.arange(8, 48, step), np.arange(8, 32, step))
bin_xy = np.column_stack([xs.ravel(), ys.ravel()])

def nearest_nucleus(xy):
    d = np.linalg.norm(nuclei[:, None, :] - xy[None, :, :], axis=2)   # (n_cells, n_bins)
    owner = d.argmin(axis=0)
    owner[d.min(axis=0) > cell_radius_um] = -1                        # outside every territory = intercellular -> label 0/drop
    return owner

owner = nearest_nucleus(bin_xy)

# each 2um bin draws a few counts (sparse) from its owner's profile; intercellular bins get ambient noise only
bin_counts = np.zeros((len(bin_xy), len(genes)))
for i, who in enumerate(owner):
    rate = cell_profiles[who] * 0.15 if who >= 0 else np.full(len(genes), 0.05)   # 0.15 keeps per-bin UMIs single-digit
    bin_counts[i] = rng.poisson(rate)

per_bin_umi = bin_counts.sum(axis=1)
print(f'2um bins: {len(bin_xy)}, median UMI/bin: {np.median(per_bin_umi):.1f} -> too sparse to cluster directly')

# RECONSTRUCT: sum bins per nucleus label (this is what b2c.bin_to_cell does)
keep = owner >= 0
cells = pd.DataFrame(bin_counts[keep], columns=genes).groupby(owner[keep]).sum()
cells['bin_count'] = pd.Series(owner[keep]).value_counts().sort_index().values   # bins absorbed per cell = QC handle
print('\nreconstructed cells (bins summed per nucleus):')
print(cells)
called = cells[genes].idxmax(axis=1)
print(f'cell identities recovered: {called.tolist()} (Epcam=epithelial, Cd3d=T cell)')

# COARSE-BIN TRAP: one 16um grid square swallows BOTH nuclei -> a mixture, not a cell
grid_um = 16.0
gx = np.floor(bin_xy[:, 0] / grid_um).astype(int)
gy = np.floor(bin_xy[:, 1] / grid_um).astype(int)
grid_id = pd.Series([f'{a}_{b}' for a, b in zip(gx, gy)])
counts_per_cell_in_grid = grid_id.groupby([grid_id, pd.Series(owner)]).size().unstack(fill_value=0)
mixed = ((counts_per_cell_in_grid[[0, 1]] > 0).sum(axis=1) >= 2).sum() if {0, 1}.issubset(counts_per_cell_in_grid.columns) else 0
print(f'\n16um grid squares mixing both cells: {mixed} -> coarse binning re-creates the deconvolution mixture')
