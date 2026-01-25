---
name: bio-hi-c-analysis-compartment-analysis
description: Detect A/B compartments from Hi-C data using cooltools and eigenvector decomposition. Identify active (A) and inactive (B) chromatin compartments from contact matrices. Use when identifying A/B compartments from Hi-C data.
tool_type: python
primary_tool: cooltools
---

# Compartment Analysis

Detect A/B compartments from Hi-C contact matrices.

## Required Imports

```python
import cooler
import cooltools
import cooltools.lib.plotting
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import bioframe
```

## Compute Compartment Eigenvectors

```python
clr = cooler.Cooler('matrix.mcool::resolutions/100000')

# Get reference genome info
view_df = bioframe.make_viewframe(clr.chromsizes)

# Compute expected values first
expected = cooltools.expected_cis(clr, view_df=view_df, ignore_diags=2)

# Compute eigenvector decomposition (compartments)
eigenvector_track = cooltools.eigs_cis(
    clr,
    view_df=view_df,
    phasing_track=None,  # Or provide GC content track
    n_eigs=3,
)

# Results are returned as a tuple (eigenvalues, eigenvectors)
eigenvalues, eigenvectors = eigenvector_track
print(f'Eigenvalues shape: {eigenvalues.shape}')
print(eigenvectors.head())
```

## Use GC Content for Phasing

```python
# GC content helps orient A/B compartments correctly
# (A compartments typically have higher GC)

# Fetch GC content
gc_track = bioframe.frac_gc(
    bioframe.make_viewframe(clr.chromsizes),
    bioframe.load_fasta('genome.fa'),
)

# Compute eigenvectors with GC phasing
eigenvalues, eigenvectors = cooltools.eigs_cis(
    clr,
    view_df=view_df,
    phasing_track=gc_track,
    n_eigs=1,
)
```

## Extract Compartment Calls

```python
# E1 (first eigenvector) defines compartments
# Positive = A (active), Negative = B (inactive)

eigenvectors['compartment'] = np.where(eigenvectors['E1'] > 0, 'A', 'B')
print(eigenvectors[['chrom', 'start', 'end', 'E1', 'compartment']].head(20))

# Count compartments
print(eigenvectors['compartment'].value_counts())
```

## Compartment Strength (Saddle Plot)

```python
# Compute saddle plot to quantify compartmentalization strength
saddle_data = cooltools.saddle(
    clr,
    expected=expected,
    eigenvector_track=eigenvectors,
    view_df=view_df,
    n_bins=50,
    vrange=(-0.5, 0.5),
)

# saddle_data contains: (saddledata, binedges)
# saddledata is the saddle matrix
saddle_matrix = saddle_data[0]
print(f'Saddle matrix shape: {saddle_matrix.shape}')
```

## Plot Saddle

```python
fig, ax = plt.subplots(figsize=(6, 6))

# Get saddle matrix (aggregate over chromosomes)
saddle_agg = np.nanmean(saddle_data[0], axis=0)

im = ax.imshow(saddle_agg, cmap='coolwarm', vmin=-1, vmax=1)
ax.set_xlabel('E1 (compartment)')
ax.set_ylabel('E1 (compartment)')
ax.set_title('Saddle plot')
plt.colorbar(im, ax=ax, label='log2(O/E)')

# Mark A and B regions
n = saddle_agg.shape[0]
ax.axhline(n/2, color='k', linewidth=0.5)
ax.axvline(n/2, color='k', linewidth=0.5)
ax.text(n*0.25, n*0.25, 'B-B', ha='center', va='center', fontsize=12)
ax.text(n*0.75, n*0.75, 'A-A', ha='center', va='center', fontsize=12)
ax.text(n*0.25, n*0.75, 'B-A', ha='center', va='center', fontsize=12)
ax.text(n*0.75, n*0.25, 'A-B', ha='center', va='center', fontsize=12)

plt.savefig('saddle_plot.png', dpi=150)
```

## Compartment Strength Score

```python
# Compute compartment strength from saddle
def compartment_strength(saddle_matrix):
    n = saddle_matrix.shape[0]
    half = n // 2

    # AA and BB corners
    AA = np.nanmean(saddle_matrix[half:, half:])
    BB = np.nanmean(saddle_matrix[:half, :half])
    AB = np.nanmean(saddle_matrix[:half, half:])
    BA = np.nanmean(saddle_matrix[half:, :half])

    # Compartment strength = (AA + BB) / (AB + BA)
    strength = (AA + BB) / 2 - (AB + BA) / 2
    return strength

strength = compartment_strength(saddle_agg)
print(f'Compartment strength: {strength:.3f}')
```

## Plot Eigenvector Track

```python
fig, ax = plt.subplots(figsize=(15, 3))

# Plot for one chromosome
chr_data = eigenvectors[eigenvectors['chrom'] == 'chr1']

# Color by compartment
colors = ['red' if e > 0 else 'blue' for e in chr_data['E1']]
ax.bar(chr_data['start'] / 1e6, chr_data['E1'], width=0.1, color=colors)

ax.axhline(0, color='k', linewidth=0.5)
ax.set_xlabel('Position (Mb)')
ax.set_ylabel('E1 (compartment)')
ax.set_title('chr1 compartments (red=A, blue=B)')

plt.tight_layout()
plt.savefig('compartment_track.png', dpi=150)
```

## Export Compartment Calls

```python
# Save as BED file
compartment_bed = eigenvectors[['chrom', 'start', 'end', 'E1', 'compartment']].copy()
compartment_bed.to_csv('compartments.bed', sep='\t', index=False, header=False)

# Save as bedGraph
eigenvectors[['chrom', 'start', 'end', 'E1']].to_csv(
    'compartment_eigenvector.bedgraph',
    sep='\t',
    index=False,
    header=False
)
```

## Compare Compartments Between Samples

```python
# Load two samples
clr1 = cooler.Cooler('sample1.mcool::resolutions/100000')
clr2 = cooler.Cooler('sample2.mcool::resolutions/100000')

# Compute eigenvectors for both
_, eig1 = cooltools.eigs_cis(clr1, view_df=view_df, n_eigs=1)
_, eig2 = cooltools.eigs_cis(clr2, view_df=view_df, n_eigs=1)

# Merge and compare
merged = eig1.merge(eig2, on=['chrom', 'start', 'end'], suffixes=('_1', '_2'))

# Correlation
from scipy.stats import pearsonr
r, p = pearsonr(merged['E1_1'].dropna(), merged['E1_2'].dropna())
print(f'E1 correlation: r={r:.3f}, p={p:.2e}')

# Compartment switches
merged['switch'] = (merged['E1_1'] > 0) != (merged['E1_2'] > 0)
print(f'Compartment switches: {merged["switch"].sum()} bins')
```

## Correlate with Gene Expression

```python
# Load gene expression data
# Assume: gene_expr with columns ['chrom', 'start', 'end', 'expression']

# Bin genes into compartment bins
compartment_expr = eigenvectors.merge(
    gene_expr,
    on=['chrom'],
    how='left'
)
compartment_expr = compartment_expr[
    (compartment_expr['start_y'] >= compartment_expr['start_x']) &
    (compartment_expr['start_y'] < compartment_expr['end_x'])
]

# Compare expression in A vs B
a_expr = compartment_expr[compartment_expr['compartment'] == 'A']['expression']
b_expr = compartment_expr[compartment_expr['compartment'] == 'B']['expression']

print(f'A compartment expression: {a_expr.mean():.2f}')
print(f'B compartment expression: {b_expr.mean():.2f}')
```

## Related Skills

- matrix-operations - Prepare matrices for compartment analysis
- hic-visualization - Visualize compartments
- chip-seq - Correlate with histone marks
