# Metagenome Visualization - Usage Guide

## Overview

This skill covers visualization and statistical analysis of metagenomic profiles using Python (matplotlib, seaborn, scikit-learn) and R (phyloseq, vegan, ggplot2).

## Common Visualizations

| Type | Purpose |
|------|---------|
| Stacked bar | Community composition |
| Heatmap | Taxa across samples |
| PCA/PCoA | Sample clustering |
| Alpha diversity | Within-sample diversity |
| Krona chart | Interactive hierarchical |

## Installation

### Python

```bash
pip install pandas matplotlib seaborn scikit-learn scipy
```

### R

```r
BiocManager::install(c('phyloseq', 'microbiome'))
install.packages(c('vegan', 'ggplot2'))
```

### Krona

```bash
conda install -c bioconda krona
```

## Input Data Formats

### MetaPhlAn Merged Table

Tab-separated with taxonomy paths as row names, samples as columns.

### Bracken Combined Output

Tab-separated with species names and read counts.

## Quick Visualization Workflow

### Python

```python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load and process
df = pd.read_csv('merged_abundance.txt', sep='\t', index_col=0)
species = df[df.index.str.contains('s__')]
species.index = species.index.str.split('|').str[-1].str.replace('s__', '')

# Top 15 heatmap
top = species.sum(axis=1).nlargest(15).index
sns.heatmap(species.loc[top], cmap='YlOrRd')
plt.tight_layout()
plt.savefig('heatmap.png')
```

### R

```r
library(phyloseq)
library(ggplot2)

# Load into phyloseq
ps <- import_metaphlan('merged_abundance.txt')

# Bar plot
plot_bar(ps, fill = 'Species')

# PCoA
ord <- ordinate(ps, 'PCoA', 'bray')
plot_ordination(ps, ord, color = 'Group')
```

## Resources

- [phyloseq Tutorial](https://joey711.github.io/phyloseq/)
- [vegan Documentation](https://cran.r-project.org/web/packages/vegan/vegan.pdf)
