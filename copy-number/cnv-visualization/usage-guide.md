# CNV Visualization Usage Guide

## Overview

Effective CNV visualization helps identify patterns, validate calls, and communicate findings. This guide covers genome-wide profiles, region-specific plots, and cohort-level comparisons.

## Prerequisites

```bash
conda install -c bioconda cnvkit
pip install matplotlib seaborn pandas plotly
```

## Quick Start

Tell your AI agent what you want to do:
- "Create a genome-wide CNV scatter plot from my CNVkit results"
- "Generate a heatmap comparing CNV profiles across all samples"
- "Plot chromosome 17 with TP53 and BRCA1 gene annotations"
- "Create an interactive CNV plot for collaborators"

## Key Plot Types

| Plot Type | Use Case |
|-----------|----------|
| Genome-wide scatter | Overview of all CNVs |
| Chromosome zoom | Detailed view of specific regions |
| Heatmap | Compare CNVs across samples |
| Ideogram/diagram | Publication-ready summary |
| Circos | Whole-genome circular view |

## Quick CNVkit Plotting

```bash
# All-in-one plots
cnvkit.py scatter sample.cnr -s sample.cns -o sample_scatter.png
cnvkit.py diagram sample.cnr -s sample.cns -o sample_diagram.pdf

# Cohort heatmap
cnvkit.py heatmap *.cns -o cohort_heatmap.pdf
```

## Publication-Quality Plots

```python
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Style settings
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 0.5

class CNVPlotter:
    def __init__(self, cnr_file, cns_file):
        self.cnr = pd.read_csv(cnr_file, sep='\t')
        self.cns = pd.read_csv(cns_file, sep='\t')
        self.setup_genome()

    def setup_genome(self):
        '''Calculate cumulative chromosome positions.'''
        self.chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX']
        self.chrom_sizes = self.cnr.groupby('chromosome')['end'].max().to_dict()

        cumsum = 0
        self.chrom_starts = {}
        self.chrom_mids = {}
        for chrom in self.chroms:
            if chrom in self.chrom_sizes:
                self.chrom_starts[chrom] = cumsum
                self.chrom_mids[chrom] = cumsum + self.chrom_sizes[chrom] / 2
                cumsum += self.chrom_sizes[chrom]
        self.genome_size = cumsum

    def genome_plot(self, figsize=(16, 3.5)):
        '''Create genome-wide CNV plot.'''
        fig, ax = plt.subplots(figsize=figsize)

        # Add cumulative position
        self.cnr['cumpos'] = self.cnr.apply(
            lambda x: self.chrom_starts.get(x['chromosome'], 0) + x['start'], axis=1)

        # Plot points
        ax.scatter(self.cnr['cumpos'] / 1e6, self.cnr['log2'],
            s=0.5, c='#888888', alpha=0.4, rasterized=True)

        # Plot segments
        for _, seg in self.cns.iterrows():
            if seg['chromosome'] in self.chrom_starts:
                start = (self.chrom_starts[seg['chromosome']] + seg['start']) / 1e6
                end = (self.chrom_starts[seg['chromosome']] + seg['end']) / 1e6
                color = '#c0392b' if seg['log2'] > 0.2 else ('#2980b9' if seg['log2'] < -0.2 else '#27ae60')
                ax.plot([start, end], [seg['log2'], seg['log2']], c=color, lw=2, solid_capstyle='butt')

        # Chromosome shading and labels
        for i, chrom in enumerate(self.chroms):
            if chrom in self.chrom_starts:
                if i % 2:
                    ax.axvspan(self.chrom_starts[chrom]/1e6,
                        (self.chrom_starts[chrom] + self.chrom_sizes[chrom])/1e6,
                        alpha=0.05, color='black')
                ax.text(self.chrom_mids[chrom]/1e6, -2.3, chrom.replace('chr', ''),
                    ha='center', fontsize=8)

        ax.set_xlim(0, self.genome_size/1e6)
        ax.set_ylim(-2.5, 2.5)
        ax.set_xlabel('Genomic position (Mb)')
        ax.set_ylabel('Log2 copy ratio')
        ax.axhline(0, c='black', lw=0.5, ls='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        return fig, ax

    def chromosome_plot(self, chrom, genes=None, figsize=(12, 3)):
        '''Plot single chromosome with optional gene annotations.'''
        fig, ax = plt.subplots(figsize=figsize)

        cnr_chr = self.cnr[self.cnr['chromosome'] == chrom]
        cns_chr = self.cns[self.cns['chromosome'] == chrom]

        ax.scatter(cnr_chr['start']/1e6, cnr_chr['log2'], s=3, c='#888888', alpha=0.5)

        for _, seg in cns_chr.iterrows():
            color = '#c0392b' if seg['log2'] > 0.2 else ('#2980b9' if seg['log2'] < -0.2 else '#27ae60')
            ax.plot([seg['start']/1e6, seg['end']/1e6], [seg['log2'], seg['log2']], c=color, lw=3)

        if genes:
            for gene, pos in genes.items():
                ax.axvline(pos/1e6, c='orange', lw=1, ls='--', alpha=0.7)
                ax.text(pos/1e6, 2, gene, rotation=90, ha='right', fontsize=8)

        ax.set_xlabel(f'{chrom} position (Mb)')
        ax.set_ylabel('Log2 copy ratio')
        ax.set_ylim(-2.5, 2.5)
        ax.axhline(0, c='black', lw=0.5, ls='--')

        plt.tight_layout()
        return fig, ax

# Usage
plotter = CNVPlotter('sample.cnr', 'sample.cns')
fig, ax = plotter.genome_plot()
fig.savefig('cnv_profile.pdf', dpi=300)

fig, ax = plotter.chromosome_plot('chr17', genes={'TP53': 7687490, 'BRCA1': 43044295})
fig.savefig('chr17_cnv.pdf', dpi=300)
```

## Heatmap for Cohort Analysis

```python
import seaborn as sns

def cohort_heatmap(cns_files, region=None, cluster=True):
    '''Create clustered heatmap of CNVs across samples.'''
    # Load and bin data
    bin_size = 1000000  # 1Mb bins
    data = {}

    for f in cns_files:
        sample = f.split('/')[-1].replace('.cns', '')
        cns = pd.read_csv(f, sep='\t')

        # Create binned representation
        bins = []
        for _, seg in cns.iterrows():
            n_bins = int((seg['end'] - seg['start']) / bin_size) + 1
            bins.extend([seg['log2']] * n_bins)
        data[sample] = bins[:1000]  # Truncate to same length

    df = pd.DataFrame(data).T

    # Cluster samples
    if cluster:
        from scipy.cluster.hierarchy import linkage, dendrogram
        linkage_matrix = linkage(df, method='ward')

    fig, ax = plt.subplots(figsize=(14, max(4, len(cns_files) * 0.4)))
    sns.heatmap(df, cmap='RdBu_r', center=0, vmin=-1, vmax=1,
        xticklabels=False, yticklabels=True, ax=ax, cbar_kws={'label': 'Log2 ratio'})

    plt.tight_layout()
    return fig

fig = cohort_heatmap(['sample1.cns', 'sample2.cns', 'sample3.cns'])
fig.savefig('cohort_cnv_heatmap.pdf')
```

## Interactive Visualization with Plotly

```python
import plotly.graph_objects as go

def interactive_cnv_plot(cnr, cns):
    '''Create interactive CNV plot.'''
    fig = go.Figure()

    # Add bins
    fig.add_trace(go.Scattergl(
        x=cnr['start'],
        y=cnr['log2'],
        mode='markers',
        marker=dict(size=2, color='gray', opacity=0.5),
        name='Bins'
    ))

    # Add segments
    for _, seg in cns.iterrows():
        color = 'red' if seg['log2'] > 0.3 else ('blue' if seg['log2'] < -0.3 else 'green')
        fig.add_trace(go.Scatter(
            x=[seg['start'], seg['end']],
            y=[seg['log2'], seg['log2']],
            mode='lines',
            line=dict(color=color, width=3),
            name=f"{seg['chromosome']}:{seg['start']}-{seg['end']}"
        ))

    fig.update_layout(
        xaxis_title='Position',
        yaxis_title='Log2 Copy Ratio',
        yaxis_range=[-2, 2],
        showlegend=False
    )

    return fig

# Save as HTML for sharing
fig = interactive_cnv_plot(cnr, cns)
fig.write_html('cnv_interactive.html')
```

## Example Prompts

> "Create a genome-wide CNV scatter plot from my CNVkit results"

> "Generate a heatmap comparing CNV profiles across all samples in my cohort"

> "Plot chromosome 17 with TP53 and BRCA1 gene annotations highlighted"

> "Create an interactive CNV plot I can share with collaborators"

## Tips for Publication Figures

1. **Resolution**: Use 300 DPI minimum
2. **Format**: PDF for vector graphics, PNG for presentations
3. **Colors**: Use colorblind-friendly palettes
4. **Labels**: Include chromosome boundaries clearly
5. **Scale**: Show Mb or kb consistently
6. **Legend**: Include log2 ratio scale interpretation
