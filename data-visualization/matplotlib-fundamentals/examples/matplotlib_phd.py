"""Reference: matplotlib 3.8+, seaborn 0.13+ | Verify API if version differs

PhD-level matplotlib pattern encoding the four correctness traps:
(1) Type-42 font embedding, (2) OO Figure/Axes API, (3) constrained_layout,
(4) rasterized point layer for vector axes + raster cells.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Illustrative stand-ins so the recipes below run as-is; replace with real results in practice.
rng = np.random.default_rng(0)
df = pd.DataFrame({'pc1': rng.normal(size=500), 'pc2': rng.normal(size=500),
                   'cluster': rng.integers(0, 4, 500), 'log_fc': rng.normal(size=500),
                   'neg_log_p': rng.exponential(2, size=500),
                   'significance': rng.choice(['up', 'down', 'ns'], 500)})
panel_data = {k: {'x': rng.normal(size=200), 'y': rng.normal(size=200)} for k in list('abcdef')}
matrix = rng.normal(size=(30, 12))

# 1. RCPARAMS FOR PUBLICATION COMPLIANCE
mpl.rcParams.update({
    'pdf.fonttype': 42,                          # TrueType -- searchable PDFs
    'ps.fonttype': 42,                           # TrueType EPS
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,                              # Nature 5-7 pt body
    'axes.labelsize': 7,
    'axes.titlesize': 8,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'lines.linewidth': 1.0,
})

# 2. NATURE SINGLE-COLUMN SCATTER (89mm wide)
fig, ax = plt.subplots(figsize=(89 / 25.4, 70 / 25.4),
                        constrained_layout=True)
ax.scatter(df['pc1'], df['pc2'],
            c=df['cluster'].astype('category').cat.codes,
            cmap='tab10', s=8, alpha=0.7,
            edgecolors='none', rasterized=True)        # CRITICAL: raster for N>1000
ax.set_xlabel('PC1 (45.2%)')                            # variance-labeled
ax.set_ylabel('PC2 (12.1%)')
ax.spines[['top', 'right']].set_visible(False)
fig.savefig('pca.pdf')                                  # Type-42 fonts; rasterized cells

# 3. MULTI-PANEL 2x3 GRID
fig, axes = plt.subplots(2, 3, figsize=(180 / 25.4, 110 / 25.4),
                          constrained_layout=True)
panel_labels = list('abcdef')
for ax, label, (key, panel) in zip(axes.flat, panel_labels, panel_data.items()):
    ax.scatter(panel['x'], panel['y'], s=4, alpha=0.6, rasterized=True)
    ax.set_title(key, fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.text(-0.15, 1.05, label, transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='top')   # Nature panel label convention

fig.savefig('multipanel.pdf')

# 4. HEATMAP with Crameri batlow + raster
from cmcrameri import cm as cmc
fig, ax = plt.subplots(figsize=(89 / 25.4, 90 / 25.4),
                        constrained_layout=True)
vmax = np.quantile(np.abs(matrix), 0.99)               # robust 99th percentile bound
im = ax.imshow(matrix, cmap=cmc.vik, vmin=-vmax, vmax=vmax,
                aspect='auto', rasterized=True)
ax.set_xlabel('Sample')
ax.set_ylabel('Gene')
cbar = fig.colorbar(im, ax=ax, shrink=0.6, aspect=20, label='Z-score')
fig.savefig('heatmap.pdf')

# 5. SEABORN INTEGRATION -- shares matplotlib Figure/Axes
import seaborn as sns
fig, ax = plt.subplots(figsize=(89 / 25.4, 70 / 25.4),
                        constrained_layout=True)
okabe_ito = ['#E69F00', '#56B4E9', '#009E73', '#F0E442',
              '#0072B2', '#D55E00', '#CC79A7', '#000000']
sns.scatterplot(data=df, x='log_fc', y='neg_log_p',
                 hue='significance', palette=okabe_ito[:3],
                 s=8, alpha=0.7, ax=ax, rasterized=True,
                 legend='brief')
ax.spines[['top', 'right']].set_visible(False)
fig.savefig('volcano_sns.pdf')

# 6. FONT EMBED VERIFICATION
# Run after savefig:
#   pdffonts pca.pdf
# Expected: "TrueType" or "Type 1" rows; NO "Type 3" rows.

# 7. SEABORN OBJECTS GRAMMAR (ggplot-like; matplotlib 3.7+, seaborn 0.13+)
import seaborn.objects as so
plot = (so.Plot(df, x='log_fc', y='neg_log_p')
          .add(so.Dots(pointsize=2), color='significance')
          .scale(color=okabe_ito[:3])
          .label(x='log2 fold change', y='-log10(p)'))
plot.save('volcano_so.pdf')
