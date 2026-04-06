---
name: bio-pathway-kegg-pathways
description: KEGG pathway and module enrichment analysis using clusterProfiler enrichKEGG and enrichMKEGG. Use when identifying metabolic and signaling pathways over-represented in a gene list. Supports 4000+ organisms via KEGG online database.
tool_type: r
primary_tool: clusterProfiler
---

## Version Compatibility

Reference examples tested with: R stats (base), clusterProfiler 4.10+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# KEGG Pathway Enrichment

## Core Pattern

**Goal:** Identify KEGG metabolic and signaling pathways over-represented in a gene list.

**Approach:** Test for enrichment using the hypergeometric test via clusterProfiler enrichKEGG against the KEGG online database.

**"Find enriched KEGG pathways in my gene list"** → Test whether KEGG pathway gene sets are over-represented among significant genes.

```r
library(clusterProfiler)

kk <- enrichKEGG(
    gene = gene_list,           # Character vector of gene IDs
    organism = 'hsa',           # KEGG organism code
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH'
)
```

## Prepare Gene List

**Goal:** Extract significant Entrez gene IDs from DE results in the format required by enrichKEGG.

**Approach:** Filter by significance thresholds and convert gene symbols to Entrez IDs (KEGG requires NCBI Entrez).

```r
library(org.Hs.eg.db)

de_results <- read.csv('de_results.csv')
sig_genes <- de_results$gene_id[de_results$padj < 0.05 & abs(de_results$log2FoldChange) > 1]

# KEGG requires NCBI Entrez gene IDs (kegg, ncbi-geneid)
gene_ids <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gene_list <- gene_ids$ENTREZID
```

## KEGG ID Conversion

**Goal:** Convert between KEGG-specific identifiers and other gene ID formats.

**Approach:** Use bitr_kegg to map between kegg, ncbi-geneid, ncbi-proteinid, and uniprot ID types.

```r
# Convert between KEGG and other IDs
kegg_ids <- bitr_kegg(gene_list, fromType = 'ncbi-geneid', toType = 'kegg', organism = 'hsa')

# Available types: kegg, ncbi-geneid, ncbi-proteinid, uniprot
```

## Run KEGG Pathway Enrichment

**Goal:** Perform KEGG pathway over-representation analysis with customizable parameters.

**Approach:** Run enrichKEGG with specified organism, ID type, and statistical thresholds.

```r
kk <- enrichKEGG(
    gene = gene_list,
    organism = 'hsa',
    keyType = 'ncbi-geneid',    # or 'kegg'
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH',
    minGSSize = 10,
    maxGSSize = 500
)

# View results
head(kk)
results <- as.data.frame(kk)
```

## Make Results Readable

```r
# enrichKEGG does NOT have readable parameter - use setReadable
library(org.Hs.eg.db)
kk_readable <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
```

## KEGG Module Enrichment

**Goal:** Test for enrichment of KEGG modules (smaller functional units than pathways).

**Approach:** Use enrichMKEGG which tests against KEGG module definitions rather than full pathways.

```r
# KEGG modules are smaller functional units than pathways
mkk <- enrichMKEGG(
    gene = gene_list,
    organism = 'hsa',
    pvalueCutoff = 0.05
)
```

## Common Organism Codes

| Code | Organism | Notes |
|------|----------|-------|
| hsa | Human (Homo sapiens) | |
| mmu | Mouse (Mus musculus) | |
| rno | Rat (Rattus norvegicus) | |
| dre | Zebrafish (Danio rerio) | |
| dme | Fruit fly (Drosophila) | |
| cel | Worm (C. elegans) | |
| sce | Yeast (S. cerevisiae) | |
| ath | Arabidopsis thaliana | |
| eco | E. coli K-12 | Bacterial |
| pae | P. aeruginosa PAO1 | Bacterial |
| bsu | B. subtilis 168 | Bacterial |
| sau | S. aureus N315 | Bacterial |
| mtc | M. tuberculosis H37Rv | Bacterial |
| ko | KEGG Orthology | Cross-species, use with KO IDs |

KEGG covers 8,000+ organisms. Always verify the code for the specific strain:
```r
search_kegg_organism('Pseudomonas', by = 'scientific_name')
search_kegg_organism('aeruginosa', by = 'scientific_name')
```

## Background Universe (Critical)

**Goal:** Restrict KEGG enrichment to genes actually measured in the experiment.

**Approach:** Convert all tested genes to Entrez IDs and pass as the universe parameter.

Without specifying the universe, enrichKEGG uses all KEGG-annotated genes as background. This inflates significance for tissue-specific pathways (e.g., liver-expressed pathways in a liver RNA-seq experiment will appear enriched simply because liver genes are expressed and brain genes are not).

```r
# Background = all tested genes (non-NA pvalue from DE analysis)
all_tested <- de_results$gene_id[!is.na(de_results$pvalue)]
universe_ids <- bitr(all_tested, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

kk <- enrichKEGG(
    gene = gene_list,
    universe = universe_ids$ENTREZID,
    organism = 'hsa',
    pvalueCutoff = 0.05
)
```

## Extract and Export Results

**Goal:** Save KEGG enrichment results to CSV and extract genes belonging to specific pathways.

**Approach:** Convert enrichment object to data frame, export, and access pathway gene sets via the geneSets slot.

```r
# Convert to data frame
results_df <- as.data.frame(kk)

# Key columns: ID (pathway), Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID, Count

# Export
write.csv(results_df, 'kegg_enrichment_results.csv', row.names = FALSE)

# Get genes in a specific pathway
pathway_genes <- kk@geneSets[['hsa04110']]  # Cell cycle
```

## Browse KEGG Pathways

**Goal:** Visualize enriched genes overlaid on KEGG pathway diagrams.

**Approach:** Use browseKEGG for interactive browser view or pathview to generate annotated pathway images.

```r
# View pathway in browser (opens KEGG website)
browseKEGG(kk, 'hsa04110')

# Download pathway image
library(pathview)
pathview(gene.data = gene_list, pathway.id = 'hsa04110', species = 'hsa')
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| gene | required | Vector of gene IDs |
| organism | hsa | KEGG organism code |
| keyType | kegg | Input ID type |
| pvalueCutoff | 0.05 | P-value threshold |
| qvalueCutoff | 0.2 | Q-value threshold |
| pAdjustMethod | BH | Adjustment method |
| universe | NULL | Background genes |
| minGSSize | 10 | Min genes per pathway |
| maxGSSize | 500 | Max genes per pathway |
| use_internal_data | FALSE | Use local KEGG data |

## Compare Multiple Gene Lists

**Goal:** Compare KEGG pathway enrichment across multiple gene lists (e.g., upregulated vs downregulated).

**Approach:** Use compareCluster with enrichKEGG to run enrichment per group and visualize with dotplot.

```r
# Compare KEGG enrichment across groups
gene_lists <- list(
    up = up_genes,
    down = down_genes
)

ck <- compareCluster(
    geneClusters = gene_lists,
    fun = 'enrichKEGG',
    organism = 'hsa'
)

dotplot(ck)
```

## Prokaryotic / Non-Model Organism KEGG

Bacteria and non-model organisms do NOT use org.*.eg.db packages or bitr(). Bacterial genes use locus tags (e.g., PA0001 for P. aeruginosa, b0001 for E. coli) that map directly as KEGG gene IDs.

```r
# Bacterial KEGG ORA -- no bitr() or OrgDb needed
# Gene IDs should be locus tags matching the KEGG genome
kegg_bac <- enrichKEGG(
    gene = sig_locus_tags,       # e.g., c('PA0001', 'PA0612', 'PA3476')
    organism = 'pae',            # P. aeruginosa PAO1
    keyType = 'kegg',            # use locus tags directly
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH'
)

# Note: setReadable() requires an OrgDb which does not exist for most bacteria
# Instead, map gene IDs manually or use KEGG gene names from the result
```

For organisms without KEGG strain-specific annotation, use KEGG Orthology (KO) with organism = 'ko'. Map genes to KO IDs via eggNOG-mapper or BlastKOALA first.

## Multi-Condition Comparison

**Goal:** Find shared and condition-specific enriched pathways across experimental conditions.

**Approach:** Run enrichKEGG per condition, then use set operations on significant pathway IDs. Do NOT compare p-values across conditions (they depend on sample size and DE gene count).

```r
# Run enrichment per condition
kk_A <- enrichKEGG(gene = sig_genes_A, organism = 'hsa', pvalueCutoff = 0.05)
kk_B <- enrichKEGG(gene = sig_genes_B, organism = 'hsa', pvalueCutoff = 0.05)

# Set operations on enriched pathway IDs
paths_A <- as.data.frame(kk_A)$ID
paths_B <- as.data.frame(kk_B)$ID
shared <- intersect(paths_A, paths_B)
only_A <- setdiff(paths_A, paths_B)
only_B <- setdiff(paths_B, paths_A)

# Or use compareCluster for side-by-side visualization
gene_clusters <- list(ConditionA = sig_genes_A, ConditionB = sig_genes_B)
ck <- compareCluster(geneClusters = gene_clusters, fun = 'enrichKEGG', organism = 'hsa')
dotplot(ck, showCategory = 10)
```

For proper multi-contrast enrichment that avoids p-value comparison pitfalls, use the mitch package (rank-MANOVA approach).

## Notes

- **No readable parameter** - use `setReadable()` with OrgDb (eukaryotes only)
- **Requires internet** - queries KEGG database online
- **use_internal_data** - set TRUE to use cached KEGG data (may be outdated)
- **Pathway IDs** - format is organism code + 5 digits (e.g., hsa04110)
- **Licensing** - KEGG data is free for academic web browsing but bulk downloads and commercial use require a license; for reproducibility-critical work, consider Reactome or WikiPathways (fully open)
- **Background universe** - always specify; default uses all KEGG-annotated genes which inflates significance

## Related Skills

- go-enrichment - Gene Ontology enrichment analysis
- gsea - GSEA using KEGG pathways (gseKEGG)
- enrichment-visualization - Visualize KEGG results
