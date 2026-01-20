# QIIME2 Workflow Usage Guide

## Overview

QIIME2 is a comprehensive microbiome analysis platform with built-in provenance tracking, a plugin architecture, and visualization tools. It provides an alternative to the DADA2/phyloseq R workflow.

## When to Use QIIME2

- **Large studies** with many samples (better parallelization)
- **Reproducibility** requirements (automatic provenance tracking)
- **Non-R users** preferring command-line workflows
- **Integration** with QIIME2 ecosystem (dozens of plugins)

## When to Use DADA2/phyloseq

- **Custom analyses** requiring R programming
- **Integration** with other Bioconductor packages
- **Interactive exploration** in RStudio
- **Existing R workflows**

## Key Concepts

### Artifacts (.qza)
Binary files containing data + metadata + provenance. Can be visualized at https://view.qiime2.org/

### Visualizations (.qzv)
Interactive HTML visualizations viewable in browser or at https://view.qiime2.org/

### Plugins
Modular components for specific tasks (dada2, deblur, feature-classifier, diversity, etc.)

## Common Workflows

### 16S V4 Region
```bash
# Standard truncation for 250bp PE reads
qiime dada2 denoise-paired \
    --p-trunc-len-f 240 \
    --p-trunc-len-r 160 \
    ...
```

### ITS Fungal
```bash
# Use ITSxpress to extract ITS region first
qiime itsxpress trim-pair-output-unmerged \
    --i-per-sample-sequences demux.qza \
    --p-region ITS2 \
    --p-taxa F \
    --o-trimmed trimmed.qza

# Then denoise without truncation
qiime dada2 denoise-paired \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    ...
```

## Choosing Sampling Depth

```bash
# Check rarefaction curves first
qiime diversity alpha-rarefaction \
    --i-table table.qza \
    --i-phylogeny rooted-tree.qza \
    --p-max-depth 50000 \
    --m-metadata-file metadata.tsv \
    --o-visualization alpha-rarefaction.qzv
```

Choose depth where curves plateau but retains most samples.

## References

- QIIME2 Documentation: https://docs.qiime2.org/
- QIIME2 Forum: https://forum.qiime2.org/
- QIIME2 Tutorials: https://docs.qiime2.org/2024.5/tutorials/
