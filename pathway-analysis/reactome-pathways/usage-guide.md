# Reactome Pathway Enrichment - Usage Guide

## Overview
ReactomePA tests a gene list (ORA via `enrichPathway`) or a ranked gene vector (GSEA via `gsePathway`) against Reactome, a curated, peer-reviewed knowledgebase whose atomic unit is the REACTION and whose pathways are nested containers of reactions. The result is therefore one signal projected onto a hierarchy, not a list of independent findings, and the honest read is the deepest curated human pathway nodes the list over-represents, deduplicated against their ancestors, against a measured background. ReactomePA reads the local `reactome.db`, so a run is reproducible offline given the Bioconductor release - unlike KEGG and WikiPathways, which query a live database.

## Prerequisites
```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('ReactomePA', 'reactome.db', 'org.Hs.eg.db', 'clusterProfiler'))
BiocManager::install('ReactomeGSA')   # only for comparative / multi-omics analysis
```

Conceptual prerequisites:
- The input is a gene LIST (ORA) or a NAMED numeric vector sorted decreasing by a per-gene statistic (GSEA).
- ReactomePA is ENTREZ-only. `enrichPathway` and `gsePathway` have no `keyType` argument, so SYMBOL or ENSEMBL ids must be converted with `bitr` first or the result is silently empty.
- ORA needs a background `universe` of the genes actually measured; omitting it defaults the background to all ~11,200 Reactome-annotated genes and inflates significance.
- Only human is curated. ReactomePA's `organism` accepts seven values (human, rat, mouse, celegans, yeast, zebrafish, fly), all non-human ones orthology-inferred from human.
- The DE list and the ranking statistic come from differential-expression; the cross-database ORA-vs-GSEA method-selection fork lives in the category README.

## Quick Start
Tell your AI agent what you want to do:
- "Run Reactome over-representation on my significant genes against my measured background"
- "Run Reactome GSEA on my ranked DESeq2 results"
- "My Reactome top hits are a parent and its child pathways - which is the real finding?"
- "Compare Reactome pathway activity between my treated and control groups"

## Example Prompts

### Over-representation
> "I have 240 significant genes as gene symbols from a DESeq2 contrast and the ~14,000 expressed genes as the background. Convert to Entrez, run Reactome over-representation, and give me the top pathways with fold enrichment, deduplicated so I am not double-counting parent and child pathways."

### GSEA
> "I have a full DESeq2 result with the test statistic for every gene. Build the ranked vector, fix the seed, and run Reactome GSEA, then show me the leading-edge genes for the top pathways."

### Hierarchy interpretation
> "My Reactome result has 'Cell Cycle Checkpoints', 'G2/M Checkpoints', and 'G1/S Transition' all near the top. Are those separate findings, and which one should I report?"

### Inspecting one pathway
> "Draw the reaction network for my top Reactome pathway colored by log2 fold change, and give me the link to open it in the Reactome Pathway Browser."

### Comparative / multi-omics
> "I have RNA-seq counts for treated vs control. Use ReactomeGSA to find which Reactome pathways differ between the groups."

## What the Agent Will Do
1. Load the DE results and extract the significant-gene list (ORA) or build the ranked named vector (GSEA).
2. Convert gene ids to Entrez with `bitr` - mandatory, since ReactomePA has no `keyType` argument.
3. Run `enrichPathway` with the measured background as `universe` and `readable=TRUE`, or `gsePathway` with a fixed seed.
4. Deduplicate the hierarchy: report the deepest significant pathway per signal and note ancestors as context.
5. Optionally draw a local reaction network with `viewPathway` (by pathway NAME) or build a PathwayBrowser URL.
6. For comparative or multi-omics questions, route to ReactomeGSA instead of ReactomePA.

## Reactome vs KEGG

| Feature | Reactome | KEGG |
|---------|----------|------|
| Atomic unit | Reaction (typed entities, PubMed-cited) | Pathway map |
| Curation | Expert-authored, externally peer-reviewed | KEGG-team curated |
| Structure | Deep event hierarchy (parent/child double-count) | Mostly flat map list |
| Granularity | Finest (reaction level); heavier multiple-testing | Coarser (pathway level) |
| Reproducibility | Local reactome.db, pinned to the Bioconductor release | Live REST API, date-dependent |
| Metabolic depth | Good | Deeper |
| License | CC0 (fully open) | Free academic; commercial license for KEGG REST |
| Organisms (in R) | 7 in ReactomePA (DB projects to ~14-20) | 8,000+ |

## Understanding Results

`enrichResult` (ORA) columns:

| Column | Description |
|--------|-------------|
| ID | Reactome stable id (R-HSA-XXXXX) |
| Description | Pathway name |
| GeneRatio | query genes in the pathway / query genes mapped to any pathway |
| BgRatio | pathway genes in the universe / universe genes mapped (denominator ~11,230 if no universe passed) |
| RichFactor | query genes in the pathway / total genes in the pathway |
| FoldEnrichment | observed / expected fraction - read effect size here, do not compute it by hand |
| zScore | standardized enrichment score |
| p.adjust | BH-adjusted p-value (report this, not raw pvalue) |
| qvalue | q-value |
| geneID | genes in the pathway (symbols when readable=TRUE) |
| Count | number of query genes in the pathway |

`gseaResult` (GSEA) adds `setSize`, `enrichmentScore`, `NES`, `rank`, `leading_edge`, and `core_enrichment`.

## Tips
- ReactomePA requires Entrez gene ids; convert from symbols or Ensembl with `bitr()` first - there is no `keyType` argument, so non-Entrez input returns nothing silently.
- Always pass a background `universe` (the genes actually measured); the implicit ~11,200-gene Reactome background otherwise inflates significance.
- Read the hierarchy, not the raw row order: a parent and its child enrich on the same genes, so report the deepest significant node and treat its ancestors as context. ReactomePA has no `simplify()` equivalent, so this is a manual call.
- `viewPathway()` takes the pathway NAME (the `Description`), NOT the R-HSA id, and it draws a LOCAL ggraph reaction network in the R graphics device - it does not open a browser. For the interactive web diagram use `browseURL('https://reactome.org/PathwayBrowser/#/<R-HSA-id>')`.
- For non-human species the annotations are orthology-inferred from human (mouse ~81% complete, less for distant species); confirm species-specific findings independently.
- The R package and the reactome.org web AnalysisService can give different p-values for the same list because of release skew and different default backgrounds; name the tool, the reactome.db version, and the background.
- For comparative (between-condition), multi-omics, or per-single-cell-cluster pathway analysis use ReactomeGSA, not ReactomePA.
- Reactome is fully open (CC0); supplement with KEGG for deeper metabolic coverage.
- See enrichment-visualization for dotplot/emapplot/cnetplot/gseaplot2 recipes.

## Related Skills

- go-enrichment - The hypergeometric test and the background-universe problem
- gsea - The GSEA running-sum engine and ranking-metric choice
- kegg-pathways - KEGG pathway/module enrichment; deeper metabolic coverage
- wikipathways - WikiPathways community-pathway enrichment (also CC0)
- enrichment-visualization - Dot/bar/cnet/emap/tree/GSEA plots of enrichment results
- differential-expression/de-results - Source of the gene list and the ranking statistic
- workflows/expression-to-pathways - End-to-end DE-to-enrichment pipeline
