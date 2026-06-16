# pathway-analysis

## Overview

Functional enrichment of gene lists and ranked gene vectors against curated gene sets (GO, KEGG, Reactome, WikiPathways, MSigDB) with R/Bioconductor: over-representation analysis (ORA), Gene Set Enrichment Analysis (GSEA), and pathway-topology analysis (SPIA). An enrichment result is a claim about a gene list versus a specific database version under a particular null, not a discovery: the BACKGROUND universe (not the gene list) decides ORA significance; GSEA needs a NAMED decreasing vector and the ranking metric IS the experiment; KEGG/WikiPathways query live databases (internet-dependent, not reproducible across releases) while GO/Reactome use local annotation; gene-ID requirements differ per method (OrgDb keyType vs kegg-id vs ENTREZ).

**Tool type:** r | **Primary tools:** clusterProfiler, ReactomePA, rWikiPathways, enrichplot

## Skills

| Skill | Description |
|-------|-------------|
| go-enrichment | GO over-representation with enrichGO; background-universe selection, ID conversion, GO-DAG redundancy reduction, RNA-seq gene-length bias (GOseq) |
| gsea | Gene Set Enrichment with gseGO/gseKEGG; named decreasing vector, ranking-metric choice, permutation type, leading edge, NES, ssGSEA/GSVA |
| kegg-pathways | KEGG pathway/module enrichment with enrichKEGG/enrichMKEGG/gseKEGG and pathway-topology analysis (SPIA/graphite); live DB and reproducibility pinning, prokaryotes, multi-condition |
| reactome-pathways | Reactome curated-pathway ORA and GSEA with ReactomePA; reaction-level detail, the nested-hierarchy double-count, ENTREZ IDs, reproducible local DB |
| wikipathways | WikiPathways enrichment with enrichWP/rWikiPathways; community curation, broad species, the live-snapshot reproducibility pin (dated GMT) |
| enrichment-visualization | enrichplot dot/bar/cnet/emap/tree/GSEA plots; the figure-is-a-modeling-choice and redundancy-collapse decision, encoding choices, required pre-steps |

## Choosing a Method

The first decision is set by the INPUT shape, not by preference (the three generations, Khatri 2012 PLoS Comput Biol 8:e1002375):

- A pre-selected gene LIST plus a defensible background (co-expression module, GWAS loci, screen hits) -> over-representation analysis (ORA, hypergeometric): go-enrichment, kegg-pathways, reactome-pathways, wikipathways. Do not binarize a ranking into a list to force ORA.
- A ranked STATISTIC for nearly all genes, no arbitrary cutoff -> Gene Set Enrichment (GSEA, running-sum): gsea. Preferred over ORA whenever a full ranking exists.
- A ranked list plus a curated signed signaling graph -> pathway topology (SPIA/graphite): kegg-pathways.

Null caveat: ORA and every clusterProfiler/fgsea preranked GSEA test a COMPETITIVE null by GENE sampling, which assumes genes are independent (they are not; co-regulated genes are positively correlated, so p-values are anti-conservative on the correlated sets of interest; Goeman and Buhlmann 2007 Bioinformatics 23:980). Treat these as a discovery screen; CAMERA (Wu and Smyth 2012 Nucleic Acids Res 40:e133) is the correlation-aware competitive test and ROAST/fry a self-contained one (both in gsea). No method is universally best (Tarca 2013 PLoS One 8:e79217; Geistlinger 2021 Brief Bioinform 22:545), so a robust finding survives a second method class.

## Method Selection

| Scenario | Method | Skill |
|----------|--------|-------|
| Ranked DE statistic for all genes, no arbitrary cutoff | GSEA | gsea |
| Pre-selected gene list (co-expression, GWAS, screens) | ORA | go-enrichment, kegg-pathways |
| Signed signaling topology + fold changes | Pathway topology (SPIA) | kegg-pathways |
| Bacterial / prokaryotic data | KEGG ORA with locus tags | kegg-pathways |
| Multiple conditions to compare | compareCluster | kegg-pathways |
| RNA-seq with gene length bias | GOseq | go-enrichment |
| Reaction-level, peer-reviewed, offline-reproducible | Reactome ORA/GSEA | reactome-pathways |
| Community-curated, broad species coverage | WikiPathways | wikipathways |

## Trustworthy-Result Checklist

Run before reporting any enrichment; each item maps to a documented failure (Wijesooriya 2022 PLoS Comput Biol 18:e1009935; Reimand 2019 Nat Protoc 14:482):

1. Correct testable-gene universe: for ORA, set `universe=` to the genes that passed the SAME filter that produced the list (e.g. the filtered DESeq2/edgeR rownames), never the genome; the default measures expression bias, not enrichment.
2. FDR-correct and report `p.adjust`/`qvalue`, never nominal `pvalue` (`pAdjustMethod='BH'` is valid under positive dependence).
3. Redundancy != replication: 40 significant GO terms are usually ~3 stories told 40 times (GO true-path rule, pathway overlap); collapse with `simplify`/EnrichmentMap and report term clusters (enrichment-visualization).
4. Leading-edge / multifunctionality check: if 3-5 hub genes explain the top ~20 sets, that is one finding wearing twenty pathway costumes (Gillis and Pavlidis 2011 PLoS One 6:e17258).
5. Version and parameter reporting: record package versions, the database/ontology release + access date (annotations drift; Tomczak 2018 Sci Rep 8:5115), the ranking metric, `pAdjustMethod`, and the universe.
6. Treat densely-annotated-set enrichment (cancer/immune/apoptosis everywhere) as the study-bias null (Haynes 2018 Sci Rep 8:1362); weight specific, less-studied hits and confirm with the leading edge.

## Gene-Set Databases

- GO (controlled function vocabulary, BP/MF/CC DAGs; broadest local coverage; heavy redundancy) and MSigDB (curated meta-collection; Hallmark H is the low-redundancy default) are the cross-cutting collections everything else samples.
- KEGG (manually drawn metabolic + signaling maps; topology-capable; live REST API), Reactome (expert-curated, peer-reviewed reactions; local reproducible reactome.db), WikiPathways (community-curated, broad species, versioned GMT) are the pathway sources, each owned by its skill.

## Example Prompts

- "Should I run ORA or GSEA on this result, and is my enrichment trustworthy?"
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
