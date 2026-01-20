# GSEA - Usage Guide

## Overview
Gene Set Enrichment Analysis (GSEA) tests whether genes in predefined sets show coordinated changes across conditions, using all genes ranked by expression change rather than just significant genes.

## Prerequisites
```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db'))

# For MSigDB gene sets:
install.packages('msigdbr')
```

## Quick Start
Tell your AI agent what you want to do:
- "Run GSEA on my differential expression results"
- "Find pathways with coordinated gene expression changes"
- "Use MSigDB Hallmark gene sets for GSEA analysis"

## Example Prompts
### Basic GSEA
> "Run GSEA on my DESeq2 results using log2FoldChange as the ranking statistic"

> "Perform GO biological process GSEA on all my genes ranked by expression change"

### MSigDB Gene Sets
> "Run GSEA using MSigDB Hallmark gene sets on my ranked gene list"

> "Use KEGG pathways from MSigDB for GSEA analysis"

### Ranking Statistics
> "Run GSEA using signed p-value as the ranking statistic instead of fold change"

> "Create a ranked gene list using the Wald statistic from DESeq2"

### Visualization
> "Show a GSEA running score plot for the top enriched pathway"

> "Create a ridge plot showing fold change distributions for enriched gene sets"

## What the Agent Will Do
1. Load DE results and create ranked gene list (named numeric vector)
2. Choose appropriate ranking statistic (log2FC, signed p-value, or Wald stat)
3. Convert gene IDs to Entrez format and sort by rank
4. Run gseGO(), gseKEGG(), or GSEA() with custom gene sets
5. Generate GSEA plots (running score, ridge plot, dotplot)

## GSEA vs Over-Representation

| Feature | Over-Representation | GSEA |
|---------|---------------------|------|
| Input | Gene list (significant only) | Ranked gene list (all genes) |
| Cutoff | Requires significance threshold | No arbitrary cutoff |
| Detection | Strong individual changes | Coordinated subtle changes |
| Functions | enrichGO, enrichKEGG | gseGO, gseKEGG |

## Choosing a Ranking Statistic

| Statistic | Formula | Best For |
|-----------|---------|----------|
| log2FC | log2FoldChange | Magnitude of change |
| Signed p | -log10(p) * sign(FC) | Both significance and direction |
| Wald | stat column from DESeq2 | Pre-computed statistic |
| t-statistic | From limma | Moderated statistics |

## Interpreting NES (Normalized Enrichment Score)
- Positive NES: Gene set genes tend to be upregulated
- Negative NES: Gene set genes tend to be downregulated
- |NES| > 1.5: Strong enrichment
- |NES| > 2.0: Very strong enrichment

## Tips
- GSEA uses ALL genes, not just significant ones - include the full ranked list
- Ensure the gene list is sorted in decreasing order before running
- Remove NAs from the ranked list before analysis
- The signed p-value statistic (-log10(p) * sign(FC)) often works best
- See enrichment-visualization skill for gseaplot2(), ridgeplot(), and dotplot()
- If no enriched terms, try a different ranking statistic or increase pvalueCutoff
