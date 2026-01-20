'''Compare Hi-C contact matrices between conditions'''

import cooler
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

clr1 = cooler.Cooler('condition1.mcool::resolutions/10000')
clr2 = cooler.Cooler('condition2.mcool::resolutions/10000')

print(f'Condition 1: {clr1.info["sum"]:,} contacts')
print(f'Condition 2: {clr2.info["sum"]:,} contacts')

region = 'chr1:50000000-60000000'
mat1 = clr1.matrix(balance=True).fetch(region)
mat2 = clr2.matrix(balance=True).fetch(region)

target = min(clr1.info['sum'], clr2.info['sum'])
mat1_norm = mat1 * (target / clr1.info['sum'])
mat2_norm = mat2 * (target / clr2.info['sum'])

log2fc = np.log2((mat2_norm + 1) / (mat1_norm + 1))
log2fc[np.isinf(log2fc)] = np.nan

print(f'\nLog2FC range: {np.nanmin(log2fc):.2f} to {np.nanmax(log2fc):.2f}')
print(f'Mean absolute log2FC: {np.nanmean(np.abs(log2fc)):.3f}')

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

im1 = axes[0].imshow(np.log2(mat1_norm + 1), cmap='Reds', vmin=-10, vmax=-3)
axes[0].set_title('Condition 1')
plt.colorbar(im1, ax=axes[0])

im2 = axes[1].imshow(np.log2(mat2_norm + 1), cmap='Reds', vmin=-10, vmax=-3)
axes[1].set_title('Condition 2')
plt.colorbar(im2, ax=axes[1])

norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)
im3 = axes[2].imshow(log2fc, cmap='coolwarm', norm=norm)
axes[2].set_title('Log2(Cond2/Cond1)')
plt.colorbar(im3, ax=axes[2])

plt.tight_layout()
plt.savefig('differential_hic.png', dpi=150)
print('\nSaved differential_hic.png')
