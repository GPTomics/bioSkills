---
name: bio-pathway-go-enrichment
description: Gene Ontology over-representation analysis using clusterProfiler enrichGO. Use when identifying biological functions enriched in a gene list from differential expression or other analyses. Supports all three ontologies (BP, MF, CC), multiple ID types, and customizable statistical thresholds.
tool_type: r
primary_tool: clusterProfiler
---

## Version Compatibility

Reference examples tested with: R stats (base), clusterProfiler 4.10+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# GO Over-Representation Analysis

## When to Use ORA vs GSEA

| Scenario | Method | Why |
|----------|--------|-----|
| Clear DE gene list with arbitrary cutoff (padj + FC) | ORA, but consider GSEA instead | ORA discards magnitude; GSEA uses all genes ranked by statistic |
| Genes from co-expression module, GWAS loci, screen hits | ORA | No ranking available; ORA is appropriate |
| All genes with DE statistics available | GSEA (gseGO) | Avoids arbitrary cutoff; detects subtle coordinated changes |
| Very few DE genes (< 20) | GSEA | ORA has no power with small lists |
| RNA-seq with known length bias | GOseq (goseq package) | Standard ORA ignores length bias; longer genes are more likely DE |

ORA converts continuous measures into binary (significant/not), losing information. When in doubt, run both ORA and GSEA and compare.

## Core Pattern

**Goal:** Identify enriched Gene Ontology terms in a gene list from differential expression or similar analyses.

**Approach:** Test for over-representation of GO terms using the hypergeometric test via clusterProfiler enrichGO.

**"Run GO enrichment on my gene list"** → Test whether biological process, molecular function, or cellular component terms are over-represented among significant genes.

```r
library(clusterProfiler)
library(org.Hs.eg.db)  # Human - change for other organisms

ego <- enrichGO(
    gene = gene_list,           # Character vector of gene IDs
    OrgDb = org.Hs.eg.db,       # Organism annotation database
    keyType = 'ENTREZID',       # ID type: ENSEMBL, SYMBOL, ENTREZID, etc.
    ont = 'BP',                 # BP, MF, CC, or ALL
    pAdjustMethod = 'BH',       # p-value adjustment method
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)
```

## Prepare Gene List from DE Results

**Goal:** Extract significant gene IDs from differential expression results and convert to the format required by enrichGO.

**Approach:** Filter DE results by adjusted p-value and fold change, then convert gene symbols to Entrez IDs using bitr.

```r
library(dplyr)

de_results <- read.csv('de_results.csv')

sig_genes <- de_results %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    pull(gene_id)

# If using gene symbols, convert to Entrez IDs
gene_ids <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gene_list <- gene_ids$ENTREZID
```

## ID Conversion with bitr

**Goal:** Convert between gene identifier types (Ensembl, Symbol, Entrez) for compatibility with enrichment tools.

**Approach:** Use clusterProfiler bitr to map between ID types using organism annotation databases.

```r
# Check available key types
keytypes(org.Hs.eg.db)

# Convert between ID types
converted <- bitr(genes, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

# Multiple output types
converted <- bitr(genes, fromType = 'SYMBOL', toType = c('ENTREZID', 'ENSEMBL'), OrgDb = org.Hs.eg.db)
```

## Background Universe (Critical)

**Goal:** Improve enrichment specificity by restricting the background to genes actually tested in the experiment.

**Approach:** Pass all expressed genes (not just significant ones) as the universe parameter to enrichGO.

The background must be genes that *could have* appeared in the list. Getting this wrong is the single most common ORA error (95% of published analyses fail to specify an appropriate background). Using the whole genome (~20,000 genes) when only 12,000 were expressed inflates significance for tissue-specific pathways.

| Experiment Type | Correct Background |
|----------------|-------------------|
| RNA-seq | All genes with detectable expression (e.g., > 1 CPM in >= N samples) |
| Microarray | All probes on the array (mapped to genes) |
| Proteomics | All detected proteins |
| Targeted panel | Only genes on the panel |

```r
# Background = all genes that were tested (NOT the full genome)
# For DESeq2: genes with non-NA pvalue survived independent filtering
all_tested <- de_results$gene_id[!is.na(de_results$pvalue)]
universe_ids <- bitr(all_tested, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

ego <- enrichGO(
    gene = gene_list,
    universe = universe_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = 'ENTREZID',
    ont = 'BP',
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05
)
```

**Warning:** clusterProfiler silently drops unannotated genes from the background. To prevent this: `options(enrichment_force_universe = TRUE)` before running enrichGO.

## All Three Ontologies

```r
# Run all ontologies at once
ego_all <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = 'ENTREZID',
    ont = 'ALL',  # BP, MF, and CC combined
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05
)

# Results include ONTOLOGY column
head(as.data.frame(ego_all))
```

## Make Results Readable

```r
# Convert Entrez IDs to gene symbols in results
ego_readable <- setReadable(ego, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')

# Or use readable = TRUE directly (only works with ENTREZID input)
ego <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = 'ENTREZID',
    ont = 'BP',
    readable = TRUE  # Converts to symbols
)
```

## Extract and Export Results

```r
# View top results
head(ego)

# Convert to data frame
results_df <- as.data.frame(ego)

# Key columns: ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count

# Export to CSV
write.csv(results_df, 'go_enrichment_results.csv', row.names = FALSE)

# Filter for specific criteria
sig_terms <- results_df[results_df$p.adjust < 0.01 & results_df$Count >= 5, ]
```

## Simplify Redundant Terms

**Goal:** Remove highly similar GO terms to reduce redundancy in enrichment results.

**Approach:** Cluster GO terms by semantic similarity and retain representative terms using the simplify function.

GO terms form a DAG (directed acyclic graph), not a flat list. If "mitotic cell cycle" is enriched, parent terms ("cell cycle", "cell cycle process") will also be enriched because they contain supersets of the same genes. Always simplify before interpretation.

```r
# Remove redundant GO terms (keeps representative terms)
ego_simplified <- simplify(ego, cutoff = 0.7, by = 'p.adjust', select_fun = min)

# measure options: 'Wang' (default, graph-based, stable across releases),
# 'Resnik', 'Lin', 'Jiang', 'Rel' (IC-based, depend on annotation version)
ego_simplified <- simplify(ego, cutoff = 0.7, measure = 'Wang')
```

**Limitations:** `simplify()` does NOT work with `ont='ALL'` -- run BP, MF, CC separately. Cutoff 0.7 is a reasonable default; lower retains more terms, higher is more aggressive.

## Different Organisms

```r
# Mouse
library(org.Mm.eg.db)
ego_mouse <- enrichGO(gene = genes, OrgDb = org.Mm.eg.db, ont = 'BP')

# Zebrafish
library(org.Dr.eg.db)
ego_zfish <- enrichGO(gene = genes, OrgDb = org.Dr.eg.db, ont = 'BP')

# Yeast
library(org.Sc.sgd.db)
ego_yeast <- enrichGO(gene = genes, OrgDb = org.Sc.sgd.db, ont = 'BP', keyType = 'ORF')
```

## Group GO Terms by Ancestor

**Goal:** Classify genes by broad GO slim categories for a high-level functional overview.

**Approach:** Use groupGO to assign genes to GO terms at a specific hierarchy level.

```r
# Classify genes by GO slim categories
ggo <- groupGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    ont = 'BP',
    level = 3,  # GO hierarchy level
    readable = TRUE
)
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| gene | required | Vector of gene IDs |
| OrgDb | required | Organism database |
| keyType | ENTREZID | Input ID type |
| ont | BP | BP, MF, CC, or ALL |
| pvalueCutoff | 0.05 | P-value threshold |
| qvalueCutoff | 0.2 | Q-value (FDR) threshold |
| pAdjustMethod | BH | BH, bonferroni, etc. |
| universe | NULL | Background genes |
| minGSSize | 10 | Min genes per term |
| maxGSSize | 500 | Max genes per term |
| readable | FALSE | Convert to symbols |

## Interpreting Results

Always examine effect size alongside p-values. A pathway with 500 genes can achieve p < 1e-15 with a modest 1.2x fold enrichment, while a 10-gene pathway with 4x enrichment at p = 0.01 is biologically more interesting.

- **Fold enrichment** = GeneRatio / BgRatio. Values > 2 suggest strong enrichment.
- **Count**: number of query genes in the term. Very large counts (> 50) may indicate overly broad terms.
- `minGSSize=10, maxGSSize=500` filters out uninformative extremes.

## Gene ID Mapping Pitfalls

- **Many-to-many mappings**: one Ensembl gene can map to multiple Entrez IDs. Deduplicate after `bitr()` to avoid counting genes multiple times.
- **Lost genes**: if > 15% of genes fail to convert, results may be unreliable. Always report the conversion rate.
- **Best practice**: use the same ID type throughout the pipeline. Convert at the last step if possible.

## RNA-seq Gene Length Bias

In RNA-seq, longer transcripts produce more fragments, increasing statistical power to detect DE. This systematically biases ORA toward pathways enriched in long genes (extracellular matrix, cell adhesion) and against short-gene pathways (ribosomal, mitochondrial). Standard normalization (RPKM, TMM) does NOT fix this.

For length-corrected GO enrichment, use GOseq:
```r
library(goseq)
pwf <- nullp(de_vector, 'hg38', 'ensGene', bias.data = gene_lengths)
goseq_results <- goseq(pwf, 'hg38', 'ensGene', method = 'Wallenius')
```

## Related Skills

- kegg-pathways - KEGG pathway enrichment
- gsea - Gene Set Enrichment Analysis for GO
- enrichment-visualization - Visualize enrichment results
- differential-expression/de-results - Generate input gene lists
