#!/usr/bin/env python3
"""Publication-ready figure export: editable fonts, hybrid rasterization, byte-stable output."""
# Reference: matplotlib 3.8+, numpy 1.26+ | Verify API if version differs

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# Okabe-Ito Color Universal Design palette (CVD-safe categorical)
OKABE_ITO = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']

# Journal column widths in mm; design at the exact target width, never rescale later
JOURNAL_WIDTH_MM = {'nature_single': 89, 'nature_double': 183, 'cell_single': 85, 'cell_double': 174}


def set_publication_style():
    """Publication rcParams: editable TrueType fonts, constrained layout, no tight-bbox."""
    mpl.rcParams.update({
        'pdf.fonttype': 42,            # TrueType -> editable/selectable text in PDF (default 3 is not)
        'ps.fonttype': 42,            # same for EPS/PS
        'svg.fonttype': 'none',        # SVG emits real <text> referencing the font
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 9,
        'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
        'axes.linewidth': 0.5, 'lines.linewidth': 1, 'lines.markersize': 4,
        'figure.dpi': 150, 'savefig.dpi': 600,
        'figure.constrained_layout.use': True,   # reproducible spacing; do NOT use savefig bbox='tight'
        'axes.spines.top': False, 'axes.spines.right': False,
    })


# Byte-stable PDFs: drop the embedded CreationDate so identical runs match
STABLE_PDF_METADATA = {'CreationDate': None}


def mm_to_in(mm):
    return mm / 25.4


def save_all(fig, stem):
    """Save vector PDF (byte-stable) + editable SVG + raster PNG."""
    fig.savefig(f'{stem}.pdf', metadata=STABLE_PDF_METADATA)
    fig.savefig(f'{stem}.svg')
    fig.savefig(f'{stem}.png')


def example_hybrid_scatter():
    """Dense scatter: rasterize the data layer, keep axes and text vector."""
    rng = np.random.default_rng(42)
    n = 200_000
    x, y = rng.standard_normal(n), rng.standard_normal(n)
    w = mm_to_in(JOURNAL_WIDTH_MM['nature_single'])
    fig, ax = plt.subplots(figsize=(w, w * 0.85))
    ax.scatter(x, y, s=1, alpha=0.3, edgecolors='none', rasterized=True, zorder=0)  # -> embedded raster
    ax.axhline(0, color=OKABE_ITO[5], lw=0.5, zorder=2)                              # stays vector
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    save_all(fig, 'figure_hybrid_scatter')
    plt.close(fig)


def example_multipanel():
    """Multi-panel figure with bold panel labels and CVD-safe colors."""
    w = mm_to_in(JOURNAL_WIDTH_MM['nature_double'])
    fig = plt.figure(figsize=(w, w * 0.55))
    gs = fig.add_gridspec(2, 3)
    rng = np.random.default_rng(0)
    ax_a = fig.add_subplot(gs[0, :2]); ax_a.plot(rng.standard_normal(100).cumsum(), color=OKABE_ITO[5])
    ax_b = fig.add_subplot(gs[0, 2]); ax_b.bar([1, 2, 3], [4, 5, 3], color=OKABE_ITO[1])
    ax_c = fig.add_subplot(gs[1, 0]); ax_c.hist(rng.standard_normal(500), bins=20, color=OKABE_ITO[3], edgecolor='white')
    ax_d = fig.add_subplot(gs[1, 1]); ax_d.scatter(rng.random(50), rng.random(50), color=OKABE_ITO[6], s=6)
    ax_e = fig.add_subplot(gs[1, 2]); ax_e.boxplot([rng.standard_normal(30) for _ in range(4)])
    for ax, label in zip([ax_a, ax_b, ax_c, ax_d, ax_e], 'ABCDE'):
        ax.text(-0.15, 1.1, label, transform=ax.transAxes, fontsize=10, fontweight='bold', va='top')
    fig.savefig('figure_multipanel.pdf', metadata=STABLE_PDF_METADATA)
    fig.savefig('figure_multipanel.png')
    plt.close(fig)


def example_tiff_for_raster_journal():
    """TIFF with LZW for journals that demand flattened raster (Cell, PLOS)."""
    w = mm_to_in(JOURNAL_WIDTH_MM['cell_single'])
    fig, ax = plt.subplots(figsize=(w, w * 0.8))
    ax.imshow(np.random.default_rng(1).random((50, 50)), cmap='cividis', aspect='auto')  # CVD-safe sequential
    fig.savefig('figure_heatmap.tiff', dpi=300, pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)


if __name__ == '__main__':
    set_publication_style()
    example_hybrid_scatter()
    example_multipanel()
    example_tiff_for_raster_journal()
    print('Figures exported: hybrid scatter (PDF/SVG/PNG), multipanel, TIFF-LZW heatmap')
