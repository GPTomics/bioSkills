# DMR Detection - Usage Guide

## Overview

Differentially Methylated Regions (DMRs) are contiguous genomic regions showing methylation differences between conditions. Multiple methods exist: tile-based (methylKit), smoothing-based (bsseq BSmooth), and kernel-based (DMRcate).

## When to Use This Skill

- You want to find regions (not just single CpGs) with methylation changes
- You have WGBS or RRBS data from multiple conditions
- You want to identify promoters or enhancers with altered methylation
- You need to integrate methylation with gene expression data

## Method Comparison

| Method | Pros | Cons | Best For |
|--------|------|------|----------|
| methylKit tiles | Simple, fast | Fixed windows | Quick exploration |
| BSmooth | Handles low coverage | Computationally intensive | WGBS |
| DMRcate | Array-optimized | Less flexible | 450K/EPIC arrays |
| DSS | Statistical rigor | Complex | Multi-factor designs |

## Choosing Tile Size (methylKit)

| Data Type | Recommended Size |
|-----------|------------------|
| WGBS | 500-1000 bp |
| RRBS | 100-500 bp |
| Targeted | Region size |

## DMR Filtering Guidelines

| Parameter | Lenient | Moderate | Stringent |
|-----------|---------|----------|-----------|
| qvalue | < 0.1 | < 0.05 | < 0.01 |
| meth.diff | > 10% | > 25% | > 40% |
| CpGs | >= 3 | >= 5 | >= 10 |

## Workflow

### 1. Single CpG Analysis First

Run methylKit analysis to understand data quality and overall patterns.

### 2. Choose DMR Method

- Quick exploration: methylKit tiles
- Publication: BSmooth or DSS
- Array data: DMRcate

### 3. Annotate DMRs

Map DMRs to genes, promoters, CpG islands.

### 4. Integrate with Other Data

Compare with gene expression, chromatin accessibility.

## Common Issues

### Few/No DMRs Found

- Relax q-value threshold
- Reduce methylation difference cutoff
- Use smaller tile size
- Check sample quality and experimental design

### Too Many DMRs

- Increase methylation difference threshold
- Lower q-value cutoff
- Merge adjacent DMRs

## Resources

- [methylKit Bioconductor](https://bioconductor.org/packages/methylKit/)
- [bsseq Bioconductor](https://bioconductor.org/packages/bsseq/)
- [DMRcate Bioconductor](https://bioconductor.org/packages/DMRcate/)
