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

| DE Tool | Recommended Metric | Column | Notes |
|---------|-------------------|--------|-------|
| DESeq2 | Wald statistic | `stat` | Best overall for RNA-seq; combines magnitude + variance |
| DESeq2 (shrunk) | Shrunken log2FC | `log2FoldChange` | Use apeglm/ashr, NOT normal type |
| limma/voom | Moderated t-statistic | `t` | Borrows strength across genes |
| edgeR | Signed p-value | `sign(logFC) * -log10(PValue)` | No Wald-equivalent; replace p=0 with 1e-300 |
| Any | log2FC alone | `log2FoldChange` | Use only when magnitude is all that matters; noisy for low-count genes |

## Interpreting NES (Normalized Enrichment Score)
- Positive NES: Gene set genes tend to be upregulated
- Negative NES: Gene set genes tend to be downregulated
- |NES| > 1.5: Strong enrichment
- Always check FDR first (< 0.25 per Broad, < 0.05 for publication), then use NES for prioritization
- High |NES| with non-significant FDR is meaningless
- Examine the leading edge genes (core_enrichment column) to see what drives the enrichment

## Tips
- GSEA uses ALL genes, not just significant ones. Include the full ranked list
- Ensure the gene list is sorted in decreasing order before running
- Remove NAs and Inf values from the ranked list before analysis
- Use the Wald statistic from DESeq2 as the default ranking metric; use signed p-value for edgeR
- Always deduplicate the ranked list after gene ID conversion (duplicates bias enrichment scores)
- FDR threshold: use < 0.25 (Broad recommendation) or < 0.05; GSEA is less powerful than ORA so uses a more lenient threshold
- Leading edge genes (core_enrichment column) are the most actionable result. Examine them to understand what drives enrichment
- See enrichment-visualization skill for gseaplot2(), ridgeplot(), and dotplot()
- If no enriched terms, try a different ranking statistic, increase pvalueCutoff, or check for duplicate gene IDs
