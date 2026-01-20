'''Plot Hi-C contact matrix'''

import cooler
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

clr = cooler.Cooler('matrix.mcool::resolutions/10000')
print(f'Loaded at {clr.binsize}bp resolution')

region = 'chr1:50000000-60000000'
matrix = clr.matrix(balance=True).fetch(region)
print(f'Matrix shape: {matrix.shape}')

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

im1 = axes[0].imshow(matrix, cmap='Reds', norm=LogNorm(vmin=0.001, vmax=0.1))
axes[0].set_title(f'{region} (balanced)')
plt.colorbar(im1, ax=axes[0], label='Contacts')

log_matrix = np.log2(matrix + 1e-10)
im2 = axes[1].imshow(log_matrix, cmap='Reds', vmin=-10, vmax=-3)
axes[1].set_title(f'{region} (log2)')
plt.colorbar(im2, ax=axes[1], label='log2(contacts)')

plt.tight_layout()
plt.savefig('hic_matrix.png', dpi=150)
print('Saved to hic_matrix.png')
