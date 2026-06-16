# Enrichment Visualization - Usage Guide

## Overview
This skill turns an enrichment result object (an enrichResult from ORA, a gseaResult from GSEA, or a compareClusterResult) into a figure with the enrichplot package. The central decision is not which plotting function to call but how to handle gene-set REDUNDANCY: a default top-20 GO dotplot is usually one biological theme drawn twenty times, because the GO DAG and nested pathway databases guarantee that a real signal surfaces as a cluster of near-identical overlapping terms. The figure either SHOWS that redundancy as structure (emapplot/treeplot via pairwise_termsim, or EnrichmentMap) or DELETES it (simplify/REVIGO). The skill also covers keeping the NES sign for GSEA, distinguishing GeneRatio from fold enrichment, and admitting that showCategory truncates.

## Prerequisites
```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('clusterProfiler', 'enrichplot', 'org.Hs.eg.db'))
BiocManager::install(c('GOSemSim', 'ggplot2'))
install.packages(c('ggridges', 'ggarchery'))   # ridgeplot needs ggridges; goplot needs ggarchery (enrichplot Suggests-only)
```

Conceptual prerequisites and notes:
- The input is an `enrichResult` / `gseaResult` / `compareClusterResult` produced by the sibling skills (go-enrichment, gsea, kegg-pathways, reactome-pathways, wikipathways), not a raw gene list.
- `ridgeplot()` and `goplot()` depend on packages enrichplot only Suggests (`ggridges` and `ggarchery`); a stock install errors with `the package ... is required` until they are installed (above).
- enrichplot 1.25.5+ moved cnetplot/emapplot/goplot to the ggtangle backend and removed several arguments (`circular`, `colorEdge`, `cex_label_gene`, `cex_label_category`, `group_category`). Introspect with `?cnetplot` / `?emapplot` and adapt; do not copy pre-2024 tutorial arguments verbatim.
- `emapplot` and `treeplot` read the `@termsim` slot and do NOT compute it - run `pairwise_termsim()` first, every time.
- There is no `barplot` method for `gseaResult` by design; GSEA results are signed and a bar cannot carry a sign.
- The plots are offline once the objects exist; building KEGG/WikiPathways objects upstream needs internet (see those skills).

## Quick Start
Tell your AI agent what you want to do:
- "Make a dotplot of my GO enrichment, collapsing redundant terms first"
- "Show the redundant terms as clusters with an enrichment map"
- "Plot my GSEA results keeping the direction (activated vs suppressed)"
- "Show a GSEA running-score plot for my top pathway"

## Example Prompts

### Collapsing redundancy
> "My enrichGO BP result has 180 significant terms and the top-20 dotplot is full of cell-cycle synonyms. Collapse the GO-DAG redundancy with simplify(), then dotplot the top 20, and tell me how many terms survived."

> "I have an enrichResult with many overlapping terms. Run pairwise_termsim and draw an enrichment map so I can see which terms are really one biological theme, and a treeplot with five named clusters."

### Encoding and effect size
> "Make a dotplot of my GO results ordered by fold enrichment instead of GeneRatio, and explain why the two orderings differ."

> "Build a gene-concept network for my top 6 enriched terms colored by log2 fold change so I can see which genes bridge multiple terms."

### GSEA plots
> "Plot my gseaResult as a ridgeplot showing the leading-edge distribution per set, keeping direction, and a gseaplot2 for the single most significant pathway."

> "Show the running enrichment score for my top three GSEA pathways overlaid in one panel with the stats table."

### Comparison and saving
> "I ran compareCluster across up- and down-regulated genes. Make a faceted dotplot comparing the two."

> "Save my enrichment dotplot as a publication-quality PDF with a viridis color scale and a caption noting the total significant term count."

## What the Agent Will Do
1. Identify the object class (enrichResult, gseaResult, or compareClusterResult) since enrichplot dispatches on it.
2. Decide the redundancy strategy: simplify()/REVIGO to delete it, or pairwise_termsim -> emapplot/treeplot to show it; never plot raw top-20 GO without a collapse step.
3. For emapplot/treeplot, compute pairwise_termsim first (JC by default; Wang with a GOSemSimDATA object for DAG-aware GO clustering).
4. Pick the encoding deliberately: GeneRatio vs FoldEnrichment for ORA, signed NES with a diverging scale for GSEA; never barplot a gseaResult.
5. Generate the figure, chain ggplot2 modifiers as needed, and save with ggsave, reporting the total significant term count and the similarity method/min_edge settings.

## Plot-by-Class Quick Reference

| Plot | Function | Class | Owns |
|------|----------|-------|------|
| Dot plot | dotplot() | ORA + GSEA | three-channel summary; default x = GeneRatio |
| Bar plot | barplot() | ORA only | unsigned count/ratio; no GSEA method |
| Gene-concept net | cnetplot() | ORA + GSEA | shared genes across terms; direction by item color |
| Enrichment map | emapplot() | ORA + GSEA | term clusters = redundancy shown (needs pairwise_termsim) |
| Tree | treeplot() | ORA + GSEA | deterministic Ward clusters (needs pairwise_termsim) |
| Ridge | ridgeplot() | GSEA only | leading-edge metric density, direction preserved |
| Running score | gseaplot2() | GSEA only | ES curve + hit ticks + ranked metric |
| Upset | upsetplot() | ORA + GSEA | gene-overlap combinations (boxplots for GSEA) |
| GO DAG | goplot() | GO only | induced DAG subgraph |
| Heatmap | heatplot() | ORA + GSEA | gene x term matrix by fold change |

## Tips
- Start from a collapse decision, not a function. For GO ORA, simplify() before a flat dotplot, or use emapplot/treeplot to show the structure.
- For emapplot/treeplot, always run pairwise_termsim() first; use method='Wang' with a GOSemSimDATA object for GO (DAG-aware) and the default method='JC' (gene overlap) for KEGG/Reactome/custom sets.
- The default dotplot orders by GeneRatio (orderBy='x'), not p-value; the top dot is the highest GeneRatio, not the most significant. Order or color by p.adjust if significance is the message.
- GeneRatio (k/n) is not fold enrichment ((k/n)/(M/N)); use x='FoldEnrichment' when specificity is the point.
- For GSEA keep the sign: use a diverging color-by-NES dotplot, ridgeplot, or gseaplot2; never coerce a gseaResult into a bar of |NES|.
- showCategory truncates; report the total number of significant terms in the caption so the figure is not read as complete.
- Lower emapplot min_edge to diagnose redundancy: if every node connects to every node, the result is one theme.
- All enrichplot functions return ggplot objects, so chain themes/scales and use ggsave; defer generic ggplot grammar to data-visualization/ggplot2-fundamentals.
- Make the object readable (setReadable) before plotting so gene labels are symbols, not Entrez IDs.

## Related Skills

- go-enrichment - Produces the enrichResult; owns simplify()
- gsea - Produces the gseaResult; owns the enrichment score and leading edge
- kegg-pathways - KEGG results to plot
- reactome-pathways - Reactome results to plot
- wikipathways - WikiPathways results to plot
- data-visualization/ggplot2-fundamentals - Generic ggplot2 grammar for the returned objects
- workflows/expression-to-pathways - End-to-end DE-to-enrichment-to-figure pipeline
