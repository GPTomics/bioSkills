# Expression to Pathways - Usage Guide

## Overview

This workflow converts differential expression results into biological insights through functional enrichment analysis using GO, KEGG, Reactome, and GSEA.

## Prerequisites

```r
BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db', 'ReactomePA', 'enrichplot'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Run pathway enrichment on my DE genes"
- "What biological processes are enriched in my significant genes?"
- "Perform GSEA on my RNA-seq results"

## Example Prompts

### Enrichment analysis
> "Run GO enrichment on my upregulated genes"

> "Find KEGG pathways for my DE gene list"

> "What Reactome pathways are enriched?"

### GSEA
> "Run GSEA using my DESeq2 results"

> "Perform gene set enrichment analysis on ranked genes"

### Visualization
> "Create a dot plot of enriched GO terms"

> "Make an enrichment map network"

> "Show me the GSEA plot for the top pathway"

## Input Requirements

| Input | Format | Required |
|-------|--------|----------|
| Gene list | Character vector | For ORA |
| Ranked genes | Named numeric vector | For GSEA |
| Organism | OrgDb package | For ID mapping |

## What the Workflow Does

1. **ID Conversion** - Convert gene symbols to Entrez IDs
2. **GO Enrichment** - Test BP, MF, CC ontologies
3. **KEGG** - Pathway over-representation
4. **Reactome** - Curated pathway enrichment
5. **GSEA** - Ranked gene set analysis
6. **Visualization** - Dot plots, networks, bar plots

## ORA vs GSEA

| Feature | ORA | GSEA |
|---------|-----|------|
| Input | Gene list | All genes ranked |
| Cutoff | Uses DE threshold | No threshold needed |
| Information | Binary (in/out) | Uses magnitude |
| Best for | Clear DE signature | Subtle changes |

## Tips

- **Gene counts**: need 50-500 genes for reliable ORA; fewer than 20 suggests switching to GSEA
- **GSEA ranking**: use Wald statistic (DESeq2), moderated t (limma), or signed p-value (edgeR), not raw p-value alone
- **Background universe**: always specify all tested genes as universe in ORA (not the full genome). This is the #1 ORA error
- **Simplify**: use simplify() to reduce redundant GO terms (does not work with ont='ALL')
- **Multiple testing**: check adjusted p-values (padj), not raw; FDR < 0.25 is acceptable for GSEA
- **Effect sizes**: examine fold enrichment and leading edge genes, not just p-values
- **Multi-condition**: use compareCluster or mitch. Never compare p-values from separate enrichment runs
- **Bacteria**: use locus tags with keyType='kegg' directly; no bitr() or OrgDb needed
- **Deduplicate**: after gene ID conversion, remove duplicates to avoid biased enrichment scores
