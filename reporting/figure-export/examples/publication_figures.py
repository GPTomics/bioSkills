#!/usr/bin/env python3
"""
Publication-ready figure export examples
Demonstrates proper sizing, resolution, and formatting for journal submission
"""

import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Publication Default Settings
# =============================================================================

def set_publication_style():
    """Configure matplotlib for publication-quality figures."""
    plt.rcParams.update({
        # Font settings
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8,
        'axes.labelsize': 8,
        'axes.titlesize': 9,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 7,

        # Line settings
        'axes.linewidth': 0.5,
        'lines.linewidth': 1,
        'lines.markersize': 4,

        # Figure settings
        'figure.dpi': 150,  # Screen display
        'savefig.dpi': 300,  # File output
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,

        # Remove top/right spines
        'axes.spines.top': False,
        'axes.spines.right': False,
    })

set_publication_style()

# =============================================================================
# Journal Figure Sizes
# =============================================================================

# Common journal column widths (inches)
JOURNAL_WIDTHS = {
    'nature_single': 89 / 25.4,      # 89mm = 3.5"
    'nature_double': 183 / 25.4,     # 183mm = 7.2"
    'cell_single': 85 / 25.4,        # 85mm
    'cell_double': 174 / 25.4,       # 174mm
    'pnas_single': 3.42,             # 8.7cm
    'pnas_double': 7.0,              # 17.8cm
    'default_single': 3.5,
    'default_double': 7.0,
}


def create_figure(width_type='default_single', height_ratio=0.75):
    """Create figure with journal-appropriate dimensions.

    Args:
        width_type: Key from JOURNAL_WIDTHS
        height_ratio: Height as fraction of width
    """
    width = JOURNAL_WIDTHS.get(width_type, 3.5)
    height = width * height_ratio
    return plt.subplots(figsize=(width, height))


# =============================================================================
# Example: Single Panel Figure
# =============================================================================

def example_single_panel():
    """Create a single-panel publication figure."""
    fig, ax = create_figure('nature_single', height_ratio=0.8)

    # Sample data (volcano plot style)
    np.random.seed(42)
    x = np.random.randn(1000)
    y = np.random.randn(1000)
    colors = np.where((np.abs(x) > 1.5) & (np.abs(y) > 1.5), 'red', 'gray')

    ax.scatter(x, y, c=colors, s=10, alpha=0.6, edgecolors='none')
    ax.axhline(1.5, ls='--', color='gray', lw=0.5)
    ax.axhline(-1.5, ls='--', color='gray', lw=0.5)
    ax.axvline(1.5, ls='--', color='gray', lw=0.5)
    ax.axvline(-1.5, ls='--', color='gray', lw=0.5)

    ax.set_xlabel('log$_2$ Fold Change')
    ax.set_ylabel('-log$_{10}$ p-value')
    ax.set_title('Differential Expression')

    # Export multiple formats
    fig.savefig('figure_volcano.pdf')   # Vector (preferred)
    fig.savefig('figure_volcano.png')   # Raster
    fig.savefig('figure_volcano.svg')   # Editable vector

    plt.close(fig)
    print('Single panel figure saved')


# =============================================================================
# Example: Multi-panel Figure
# =============================================================================

def example_multipanel():
    """Create a multi-panel figure with labels."""
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(JOURNAL_WIDTHS['default_double'], 4))
    gs = GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.4)

    # Panel A: Large plot spanning 2 columns
    ax_a = fig.add_subplot(gs[0, :2])
    ax_a.plot(np.random.randn(100).cumsum())
    ax_a.set_xlabel('Time')
    ax_a.set_ylabel('Value')

    # Panel B: Small plot
    ax_b = fig.add_subplot(gs[0, 2])
    ax_b.bar([1, 2, 3], [4, 5, 3])

    # Panel C: Bottom row
    ax_c = fig.add_subplot(gs[1, 0])
    ax_c.hist(np.random.randn(500), bins=20, edgecolor='white')

    # Panel D
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.scatter(np.random.rand(50), np.random.rand(50))

    # Panel E
    ax_e = fig.add_subplot(gs[1, 2])
    ax_e.boxplot([np.random.randn(30) for _ in range(4)])

    # Add panel labels
    panels = [(ax_a, 'A'), (ax_b, 'B'), (ax_c, 'C'), (ax_d, 'D'), (ax_e, 'E')]
    for ax, label in panels:
        ax.text(-0.15, 1.1, label, transform=ax.transAxes,
                fontsize=10, fontweight='bold', va='top')

    fig.savefig('figure_multipanel.pdf')
    fig.savefig('figure_multipanel.png')
    plt.close(fig)
    print('Multi-panel figure saved')


# =============================================================================
# Colorblind-safe Palettes
# =============================================================================

# Recommended colorblind-safe palettes
COLORBLIND_PALETTES = {
    'categorical': ['#0077BB', '#33BBEE', '#009988', '#EE7733', '#CC3311', '#EE3377'],
    'diverging': plt.cm.RdBu,
    'sequential': plt.cm.viridis,
}


# =============================================================================
# Main
# =============================================================================

if __name__ == '__main__':
    example_single_panel()
    example_multipanel()
    print('\nAll figures exported successfully')
    print('Formats: PDF (vector), PNG (300 DPI), SVG (editable)')
