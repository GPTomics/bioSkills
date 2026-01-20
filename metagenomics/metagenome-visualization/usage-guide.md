# Metagenome Visualization - Usage Guide

## Overview
Visualize and statistically analyze metagenomic profiles using Python (matplotlib, seaborn, scikit-learn) or R (phyloseq, vegan, ggplot2).

## Prerequisites
```bash
# Python
pip install pandas matplotlib seaborn scikit-learn scipy

# R
Rscript -e "BiocManager::install(c('phyloseq', 'microbiome'))"
Rscript -e "install.packages(c('vegan', 'ggplot2'))"

# Krona for interactive charts
conda install -c bioconda krona
```

## Quick Start
Tell your AI agent what you want to do:
- "Create a stacked bar chart of community composition"
- "Make a heatmap of the top 20 species across samples"
- "Plot PCoA to visualize sample clustering"

## Example Prompts
### Composition Plots
> "Create a stacked bar plot showing phylum-level composition for all samples"

> "Make a heatmap of the top 15 most abundant species across my samples"

### Ordination and Clustering
> "Run PCoA with Bray-Curtis distance and color by treatment group"

> "Cluster my samples based on species profiles and show a dendrogram"

### Diversity Analysis
> "Calculate and plot alpha diversity (Shannon, Simpson) for each group"

> "Create rarefaction curves for all my samples"

### Interactive Visualization
> "Generate a Krona chart from my Kraken2 output"

> "Make an interactive plot where I can hover to see species names"

## What the Agent Will Do
1. Load and parse abundance data (MetaPhlAn, Bracken, or other formats)
2. Filter to relevant taxonomic levels and aggregate if needed
3. Create publication-quality visualizations
4. Save figures in requested format (PNG, PDF, SVG)

## Tips
- MetaPhlAn outputs relative abundance (sums to 100%)
- Bracken outputs read counts (normalize before comparing)
- Use log transformation for highly skewed data
- phyloseq objects integrate abundance, taxonomy, and metadata
- vegan's `vegdist` supports many distance metrics (Bray-Curtis, Jaccard, etc.)

## Common Visualizations

| Type | Purpose |
|------|---------|
| Stacked bar | Community composition |
| Heatmap | Taxa across samples |
| PCA/PCoA | Sample clustering |
| Alpha diversity | Within-sample diversity |
| Krona chart | Interactive hierarchical |

## Resources
- [phyloseq Tutorial](https://joey711.github.io/phyloseq/)
- [vegan Documentation](https://cran.r-project.org/web/packages/vegan/vegan.pdf)
