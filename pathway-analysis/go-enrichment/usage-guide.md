# GO Enrichment - Usage Guide

## Overview
Gene Ontology over-representation analysis tests whether specific GO terms appear more frequently in your gene list than expected by chance.

## Prerequisites
```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Run GO enrichment on my differentially expressed genes"
- "Find biological processes enriched in my upregulated genes"
- "Do GO analysis comparing treatment vs control"

## Example Prompts
### Basic Enrichment
> "Run GO enrichment on the significant genes in deseq2_results.csv using padj < 0.05 and log2FC > 1 cutoffs"

> "Perform GO biological process enrichment on my gene list and show the top 20 terms"

### Separate Up/Down Analysis
> "Run GO enrichment separately for upregulated and downregulated genes"

> "Compare GO terms between genes that go up vs down in treatment"

### Advanced Options
> "Run GO enrichment using all expressed genes as background universe"

> "Simplify redundant GO terms after enrichment analysis"

> "Run GO enrichment for all three ontologies (BP, MF, CC) and compare results"

## What the Agent Will Do
1. Load DE results and extract significant genes based on specified cutoffs
2. Convert gene IDs to Entrez format using bitr()
3. Run enrichGO() with appropriate ontology and parameters
4. Optionally simplify redundant terms with simplify()
5. Generate results table and visualizations

## Understanding Results

| Column | Description |
|--------|-------------|
| ID | GO term ID (GO:XXXXXXX) |
| Description | GO term name |
| GeneRatio | Genes in term / Total query genes |
| BgRatio | Background genes in term / Total background |
| pvalue | Raw p-value (hypergeometric test) |
| p.adjust | Adjusted p-value (FDR) |
| qvalue | Q-value |
| geneID | Genes in the term |
| Count | Number of genes |

## Three Ontologies

| Ontology | Code | Description |
|----------|------|-------------|
| Biological Process | BP | What the genes do |
| Molecular Function | MF | Biochemical activity |
| Cellular Component | CC | Where in cell |

## Tips
- Always specify a background universe (all tested genes, not the full genome). This is the single most common ORA error
- Use `options(enrichment_force_universe = TRUE)` to prevent clusterProfiler from dropping unannotated genes from the background
- Use simplify() to reduce redundant hierarchical terms (cutoff = 0.7 is a good start); simplify() does not work with ont='ALL', so run BP, MF, CC separately
- Run separate enrichment on up- and down-regulated genes to see directional patterns
- Examine fold enrichment (GeneRatio / BgRatio), not just p-values, to assess biological significance
- For RNA-seq data, consider GOseq (goseq package) which corrects for gene length bias
- If no terms are found, check gene ID conversion success rate and loosen thresholds
- After bitr(), deduplicate to avoid counting genes with many-to-many ID mappings multiple times
- See enrichment-visualization skill for plotting options (dotplot, cnetplot, emapplot)
- Treat enrichment results as hypothesis generation, not validation. Enrichments derived from DE results cannot "validate" those same DE results
