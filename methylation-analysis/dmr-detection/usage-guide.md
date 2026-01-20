# DMR Detection - Usage Guide

## Overview
Differentially Methylated Region (DMR) detection identifies contiguous genomic regions with methylation differences between conditions, using methods like tile-based (methylKit), smoothing-based (bsseq BSmooth), or kernel-based (DMRcate) approaches.

## Prerequisites
```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('methylKit', 'bsseq', 'DMRcate', 'DSS'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Find differentially methylated regions between my treatment groups"
- "Detect DMRs using BSmooth smoothing for my WGBS data"
- "Identify promoter regions with altered methylation"

## Example Prompts
### Method Selection
> "Which DMR method should I use for my WGBS data with low coverage?"

> "Compare tile-based vs smoothing-based DMR detection for my samples"

### methylKit Tiles
> "Find DMRs using 500bp tiles with methylKit"

> "Run tile-based DMR analysis with stringent cutoffs"

### BSmooth Analysis
> "Detect DMRs using BSmooth smoothing for my WGBS experiment"

> "Run bsseq DMR analysis with appropriate smoothing parameters"

### DMR Annotation
> "Annotate my DMRs with gene names and genomic features"

> "Find DMRs overlapping promoters and CpG islands"

### Integration
> "Compare DMRs with differentially expressed genes"

> "Identify DMRs in enhancer regions from my ChIP-seq data"

## What the Agent Will Do
1. Help select appropriate DMR method based on your data type
2. Set tile size or smoothing parameters for your analysis
3. Run DMR detection with chosen method
4. Filter DMRs by q-value, methylation difference, and number of CpGs
5. Annotate DMRs with genomic features
6. Export results for visualization or integration

## Method Comparison

| Method | Pros | Cons | Best For |
|--------|------|------|----------|
| methylKit tiles | Simple, fast | Fixed windows | Quick exploration |
| BSmooth | Handles low coverage | Computationally intensive | WGBS |
| DMRcate | Array-optimized | Less flexible | 450K/EPIC arrays |
| DSS | Statistical rigor | Complex | Multi-factor designs |

## DMR Filtering Guidelines

| Parameter | Lenient | Moderate | Stringent |
|-----------|---------|----------|-----------|
| qvalue | < 0.1 | < 0.05 | < 0.01 |
| meth.diff | > 10% | > 25% | > 40% |
| CpGs | >= 3 | >= 5 | >= 10 |

## Tips
- Run single-CpG analysis first to understand data quality before DMR detection
- Tile size recommendations: WGBS 500-1000bp, RRBS 100-500bp
- BSmooth works better with low coverage data due to smoothing
- For publication, use BSmooth or DSS rather than simple tiles
- If few DMRs found: relax thresholds, use smaller tiles, check sample quality
- If too many DMRs: increase meth.diff threshold, lower q-value cutoff
- Always annotate DMRs to promoters, CpG islands, and gene bodies for interpretation
- Consider integrating with expression data to prioritize functionally relevant DMRs
