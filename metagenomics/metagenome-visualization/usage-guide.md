# Metagenome Visualization - Usage Guide

## Overview
This skill turns a shotgun profiler table (MetaPhlAn relative abundance, Bracken counts, HUMAnN function tables) into honest figures and defensible community statistics with phyloseq, vegan, microViz, and Python. The central idea: an ordination, a stacked bar, and a diversity number are modeling choices - the distance, transform, rarefaction, filter, and top-N cutoff can each create or erase a result, so they must be declared and shown not to fabricate the conclusion. The MetaPhlAn-percent vs Bracken-counts fork decides what is valid downstream.

## Prerequisites
```bash
# R is primary for community statistics
Rscript -e "BiocManager::install(c('phyloseq', 'microbiome', 'ALDEx2', 'ANCOMBC'))"
Rscript -e "install.packages(c('vegan', 'microViz'))"
# Python for plotting/wrangling
pip install pandas matplotlib seaborn scikit-bio
conda install -c bioconda krona && ktUpdateTaxonomy.sh
```

Conceptual prerequisites:
- MetaPhlAn = percent (cannot rarefy; SGB tree so UniFrac available); Bracken = counts (rarefaction/count models valid; no tree).
- Bray-Curtis is the default but compositionally incoherent; Aitchison (CLR) is correct.
- Shotgun richness tracks database size and false positives; prefer Hill q=1/q=2.
- Pair every PERMANOVA with betadisper; use a consensus of >=2 differential-abundance tools.

## Quick Start
Tell your AI agent what you want to do:
- "Make an honest stacked bar that labels how many taxa are hidden in Other"
- "Run a compositional ordination (CLR-PCA / Aitchison) and check it against Bray-Curtis"
- "Test differential abundance with two compositional tools and report the intersect"
- "Compute Hill-number diversity and show a rarefaction curve"

## Example Prompts

### Composition and ordination
> "Build a phyloseq object from my MetaPhlAn table, run a CLR-PCA ordination colored by group, and show whether the separation survives a Bray-Curtis PCoA."

### PERMANOVA done right
> "Run adonis2 on Bray-Curtis distance by treatment, and pair it with betadisper so I know whether this is a location shift or a dispersion difference."

### Differential abundance consensus
> "Test differential abundance with ALDEx2 and ANCOM-BC on my Bracken counts, prevalence-filter first, BH-correct, and report the taxa both tools agree on."

### Diversity honesty
> "Compute Hill numbers q=0/1/2, show a rarefaction curve, and tell me how much my richness depends on the Kraken confidence threshold."

## What the Agent Will Do
1. Parse the profiler table, filter to species (or functions), and note whether it is counts or percentages.
2. Choose a transform and distance for a stated reason and show the conclusion survives an alternative.
3. Pair PERMANOVA with a dispersion test and report both.
4. Run >=2 compositional DA tools, prevalence-filter, BH-correct, and report the intersect.
5. Report diversity as Hill numbers, gate richness behind confidence filtering, and show a rarefaction curve.
6. Label what a stacked bar hides (Other count/percent, n, absolute load).

## Tips
- Krona shows the full drillable hierarchy - often more honest than a top-10 bar.
- Never run Chao1/Observed on MetaPhlAn relative abundances; richness estimators need integer counts.
- CLR before PCA (Aitchison-PCA), never StandardScaler on raw relative abundance.
- UniFrac/Faith's PD need a tree - available for MetaPhlAn SGBs, not Bracken/HUMAnN.
- Name every DA tool; single-tool hits are tentative.

## Common Visualizations

| Type | Purpose | Caveat |
|------|---------|--------|
| Krona | interactive hierarchy | the honest "what is in it" |
| Stacked bar | composition | relative race; hides absolute load |
| CLR heatmap | taxa across samples | transform is a modeling choice |
| CLR-PCA / PCoA | sample structure | declare metric; check an alternative |
| Hill numbers | within-sample diversity | richness is a database readout |

## Related Skills

- metaphlan-profiling - Generates the relative-abundance table with the SGB tree
- kraken-classification - Generates Kraken/Bracken count input
- abundance-estimation - Compositional theory, normalization, absolute load
- functional-profiling - HUMAnN function tables tested with the same DA logic
- microbiome/diversity-analysis - Amplicon/QIIME2 diversity for ASV input
- data-visualization/ggplot2-fundamentals - Generic plotting primitives
- workflows/metagenomics-pipeline - End-to-end shotgun analysis

## Resources
- [phyloseq tutorial](https://joey711.github.io/phyloseq/)
- [microViz](https://david-barnett.github.io/microViz/)
- [vegan documentation](https://cran.r-project.org/web/packages/vegan/vegan.pdf)
