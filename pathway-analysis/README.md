# pathway-analysis

## Overview

Functional enrichment of gene lists and ranked gene vectors against curated gene sets (GO, KEGG, Reactome, WikiPathways, MSigDB) with R/Bioconductor: over-representation analysis (ORA), Gene Set Enrichment Analysis (GSEA), and pathway-topology analysis (SPIA). Decision-grade framing: the first decision is the generation (ORA vs GSEA vs topology), set by the input not by preference; the BACKGROUND universe (not the gene list) decides ORA significance; GSEA needs a NAMED decreasing vector and the ranking metric IS the experiment; KEGG/WikiPathways query live databases (internet-dependent, not reproducible across releases) while GO/Reactome use local annotation; gene-ID requirements differ per method (OrgDb keyType vs kegg-id vs ENTREZ).

**Tool type:** r | **Primary tools:** clusterProfiler, ReactomePA, rWikiPathways, enrichplot

## Skills

| Skill | Description |
|-------|-------------|
| enrichment-foundations | Method-selection spine: the ORA-vs-FCS-vs-topology generation fork, competitive vs self-contained nulls, the master decision tree, and the is-my-result-trustworthy layer (background universe, FDR/version reporting, redundancy != replication, annotation bias) |
| go-enrichment | GO over-representation with enrichGO; background-universe selection, ID conversion, GO-DAG redundancy reduction, RNA-seq gene-length bias (GOseq) |
| gsea | Gene Set Enrichment with gseGO/gseKEGG; named decreasing vector, ranking-metric choice, permutation type, leading edge, NES, ssGSEA/GSVA |
| kegg-pathways | KEGG pathway/module enrichment with enrichKEGG/enrichMKEGG/gseKEGG and pathway-topology analysis (SPIA/graphite); live DB and reproducibility pinning, prokaryotes, multi-condition |
| reactome-pathways | Reactome curated-pathway ORA and GSEA with ReactomePA; reaction-level detail, the nested-hierarchy double-count, ENTREZ IDs, reproducible local DB |
| wikipathways | WikiPathways enrichment with enrichWP/rWikiPathways; community curation, broad species, the live-snapshot reproducibility pin (dated GMT) |
| enrichment-visualization | enrichplot dot/bar/cnet/emap/tree/GSEA plots; the figure-is-a-modeling-choice and redundancy-collapse decision, encoding choices, required pre-steps |

## Method Selection

| Scenario | Method | Skill |
|----------|--------|-------|
| Unsure which method at all / is my result trustworthy | The generation fork + trustworthiness checklist | enrichment-foundations |
| Ranked DE statistic for all genes, no arbitrary cutoff | GSEA | gsea |
| Pre-selected gene list (co-expression, GWAS, screens) | ORA | go-enrichment, kegg-pathways |
| Signed signaling topology + fold changes | Pathway topology (SPIA) | kegg-pathways |
| Bacterial / prokaryotic data | KEGG ORA with locus tags | kegg-pathways |
| Multiple conditions to compare | compareCluster | kegg-pathways |
| RNA-seq with gene length bias | GOseq | go-enrichment |

## Example Prompts

- "Should I run ORA or GSEA on this result, and is my enrichment trustworthy?"
- "Pick the right enrichment method for my data and check the background universe"
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
BiocManager::install(c('ReactomePA', 'rWikiPathways', 'msigdbr'))
# For gene length bias correction in RNA-seq:
BiocManager::install('goseq')
# For KEGG pathway-topology analysis (SPIA):
BiocManager::install(c('SPIA', 'graphite'))
```

## Related Skills

- **differential-expression** - Generate gene lists and statistics for enrichment
- **single-cell** - Marker genes can be analyzed with pathway enrichment
- **database-access** - Fetch gene annotations from NCBI
- **workflows** - expression-to-pathways orchestrates the full DE-to-enrichment pipeline
