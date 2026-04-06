# pathway-analysis

## Overview

Functional enrichment and pathway analysis using R/Bioconductor. Supports over-representation analysis (ORA) and Gene Set Enrichment Analysis (GSEA) across Gene Ontology, KEGG, Reactome, and WikiPathways databases. Includes guidance for prokaryotic organisms, multi-condition comparison, and common methodological pitfalls.

**Tool type:** r | **Primary tools:** clusterProfiler, ReactomePA, rWikiPathways, enrichplot

## Skills

| Skill | Description |
|-------|-------------|
| go-enrichment | Gene Ontology over-representation analysis with enrichGO; includes ORA vs GSEA decision guidance, background universe selection, gene length bias correction |
| kegg-pathways | KEGG pathway enrichment with enrichKEGG and enrichMKEGG; includes prokaryotic organism support and multi-condition comparison |
| reactome-pathways | Reactome pathway enrichment with ReactomePA; peer-reviewed curated pathways with reaction-level detail |
| wikipathways | WikiPathways enrichment with enrichWP and rWikiPathways; community-curated open-source pathways for 30+ species |
| gsea | Gene Set Enrichment Analysis with gseGO, gseKEGG; ranking metric selection, leading edge interpretation, NES caveats |
| enrichment-visualization | Dot plots, bar plots, enrichment maps, cnetplots, GSEA plots; visualization selection guidance and common mistakes |

## Method Selection

| Scenario | Method | Skill |
|----------|--------|-------|
| Ranked DE results for all genes | GSEA | gsea |
| Gene list from co-expression, GWAS, screens | ORA | go-enrichment, kegg-pathways |
| Bacterial / prokaryotic data | KEGG ORA with locus tags | kegg-pathways |
| Multiple conditions to compare | compareCluster or mitch | kegg-pathways |
| RNA-seq with gene length bias | GOseq | go-enrichment |

## Example Prompts

- "Run GO enrichment on my differentially expressed genes"
- "Find enriched biological processes for these genes"
- "What molecular functions are over-represented in my gene list?"
- "Find enriched KEGG pathways for my gene set"
- "What pathways are active in my differentially expressed genes?"
- "Run KEGG module enrichment analysis"
- "Run KEGG enrichment on my P. aeruginosa DE results"
- "Run Reactome pathway enrichment on my genes"
- "Find enriched Reactome pathways for my DEGs"
- "Run WikiPathways enrichment analysis"
- "Run GSEA on my ranked gene list"
- "Perform gene set enrichment analysis using GO terms"
- "Run GSEA with KEGG pathways"
- "Compare enriched pathways between treatment and control conditions"
- "Create a dot plot of my enrichment results"
- "Make an enrichment map showing term relationships"
- "Show a gene-concept network for top pathways"
- "Create a GSEA running score plot"

## Requirements

```r
BiocManager::install(c('clusterProfiler', 'enrichplot', 'org.Hs.eg.db'))
BiocManager::install(c('ReactomePA', 'rWikiPathways'))
# For gene length bias correction in RNA-seq:
BiocManager::install('goseq')
```

## Related Skills

- **differential-expression** - Generate gene lists and statistics for enrichment
- **single-cell** - Marker genes can be analyzed with pathway enrichment
- **database-access** - Fetch gene annotations from NCBI
- **workflows** - expression-to-pathways orchestrates the full DE-to-enrichment pipeline
