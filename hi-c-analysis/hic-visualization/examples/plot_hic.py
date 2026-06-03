'''Render a Hi-C region three ways (balanced / observed-over-expected / log2-ratio) and an APA pileup.

The colorscale is where Hi-C figures lie: balanced shows TADs but hides compartments/loops,
which only the symmetric log2(O/E) view reveals. Every panel records normalization + scale +
resolution so the figure is reproducible from its legend.
'''
# Reference: cooler 0.10+, cooltools 0.7+, bioframe 0.7+, matplotlib 3.8+ | Verify API if version differs

import cooler
import cooltools
import bioframe
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

try:
    import cooltools.lib.plotting   # registers the canonical 'fall' Hi-C colormap (needs matplotlib < 3.9 with cooltools 0.7.x)
    HIC_CMAP = 'fall'
except ImportError:
    HIC_CMAP = 'afmhot_r'           # stock white-yellow-red-black fallback when 'fall' is unavailable

RESOLUTION = 10_000        # 10kb: TAD/loop scale; use 100-500kb for compartments
VMAX_PCTL = 99.5           # off-diagonal percentile for vmax; heavy-tailed counts need clipping, not data-max
OE_CLIP = 2.0              # symmetric |log2(O/E)| limit; ~2 is the usual readable range
RATIO_CLIP = 2.0           # symmetric log2-ratio limit; vmin=-vmax keeps white at no-change
APA_FLANK = 100_000        # +/-100kb APA window; too small contaminates the corner background
PSEUDOCOUNT = 1e-5         # tames divide-by-small off-diagonal in the log2-ratio

region = ('chr1', 50_000_000, 60_000_000)

clr = cooler.Cooler(f'matrix.mcool::/resolutions/{RESOLUTION}')
view_df = bioframe.make_viewframe(clr.chromsizes)
expected = cooltools.expected_cis(clr, view_df=view_df, nproc=4)

balanced = clr.matrix(balance=True).fetch(region)
vmax = np.nanpercentile(balanced[balanced > 0], VMAX_PCTL)
exp_by_diag = expected.query('region1 == @region[0]')['balanced.avg'].to_numpy()   # expected per separation
ii, jj = np.indices(balanced.shape)
oe = balanced / exp_by_diag[np.abs(ii - jj)]   # each pixel / its distance-matched expected

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

fall = plt.get_cmap(HIC_CMAP).copy(); fall.set_bad('lightgray')   # masked bins shown, not white
im0 = axes[0].matshow(balanced, norm=LogNorm(vmin=vmax * 1e-3, vmax=vmax), cmap=fall)
axes[0].set_title(f'balanced (ICE), {RESOLUTION // 1000}kb, vmax=p{VMAX_PCTL}')
fig.colorbar(im0, ax=axes[0], fraction=0.046)

im1 = axes[1].matshow(np.log2(oe), cmap='coolwarm', vmin=-OE_CLIP, vmax=OE_CLIP)   # symmetric: 0 stays white
axes[1].set_title('log2(obs/exp) -- compartments + loops')
fig.colorbar(im1, ax=axes[1], fraction=0.046)

clr2 = cooler.Cooler(f'condition2.mcool::/resolutions/{RESOLUTION}')   # already balanced + depth-matched
balanced2 = clr2.matrix(balance=True).fetch(region)
ratio = np.log2((balanced + PSEUDOCOUNT) / (balanced2 + PSEUDOCOUNT))
im2 = axes[2].matshow(ratio, cmap='RdBu_r', vmin=-RATIO_CLIP, vmax=RATIO_CLIP)
axes[2].set_title('log2(cond1/cond2)')
fig.colorbar(im2, ax=axes[2], fraction=0.046)

fig.tight_layout()
fig.savefig('hic_transforms.png', dpi=150)

loops = pd.read_csv('loops.bedpe', sep='\t')   # chrom1,start1,end1,chrom2,start2,end2
stack = cooltools.pileup(clr, loops, view_df=view_df, expected_df=expected, flank=APA_FLANK)
apa = np.nanmean(stack, axis=0)   # stack is (n_features, D, D) -> average over the loop axis
center = apa.shape[0] // 2
apa_score = apa[center, center] / np.nanmean(apa[-3:, :3])   # center / lower-left 3x3 background (Rao 2014)

fig2, ax = plt.subplots(figsize=(5, 5))
im = ax.matshow(np.log2(apa), cmap='coolwarm', vmin=-1, vmax=1)
ax.set_title(f'APA score {apa_score:.2f}')
fig2.colorbar(im, ax=ax, fraction=0.046)
fig2.savefig('apa_pileup.png', dpi=150)
