# KEGG Pathway and Topology Enrichment - Usage Guide

## Overview
KEGG enrichment tests genes against KEGG's curated pathway and module gene sets across all three generations of pathway analysis: over-representation (enrichKEGG, enrichMKEGG), ranked GSEA (gseKEGG), and signed-topology perturbation (SPIA/graphite). KEGG is the database that owns the third generation because it ships signed directed signaling topology (KGML), letting SPIA propagate fold-changes through the wiring rather than treating a pathway as an unordered gene set. The caveat that breaks reproducibility: a KEGG result is a timestamped query against a live, partially-paywalled REST API, so it is irreproducible unless the release is pinned.

## Prerequisites
```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db'))
install.packages('gson')                              # snapshot pinning for reproducibility
BiocManager::install(c('SPIA', 'graphite'))           # signed-topology perturbation
BiocManager::install('pathview')                      # KEGG-map overlay
```

Conceptual prerequisites:
- The input is a gene LIST (ORA), a ranked named-decreasing vector (GSEA), or a named fold-change vector plus a universe (SPIA). The DE list and fold-changes come from differential-expression.
- Gene IDs must match what KEGG expects: Entrez (keyType='ncbi-geneid') for model eukaryotes, or locus tags (keyType='kegg') for prokaryotes. ENSEMBL/SYMBOL passed to enrichKEGG return zero hits. There is no OrgDb and no Entrez==KEGG identity for most bacteria.
- The background universe must be the genes that could have been called DE, not all KEGG-annotated genes.
- enrichKEGG/enrichMKEGG/gseKEGG/SPIA query the KEGG REST API at runtime: they need internet and are NOT reproducible across KEGG releases. For any reported result, pin with a gson snapshot and record the access date.
- SPIA is signaling-only; it is undefined for metabolic maps.

## Quick Start
Tell your AI agent what you want to do:
- "Run KEGG pathway enrichment on my significant genes with the measured genes as background"
- "Score signed KEGG pathway perturbation from my fold-changes with SPIA"
- "Run KEGG enrichment on my bacterial DE list using locus tags"
- "Pin the KEGG release so my enrichment is reproducible"
- "Compare KEGG enrichment between my up- and down-regulated genes"

## Example Prompts

### KEGG over-representation
> "I have 240 significant genes from a DESeq2 contrast as SYMBOLs and the ~13,000 expressed genes as background. Convert both to Entrez, run KEGG pathway ORA for human with the measured universe, and give me the top pathways by adjusted p-value with fold enrichment."

### Signed-topology perturbation
> "I have a human DE list with log2 fold-changes and a universe. Run SPIA so direction and network position are used, and tell me which signaling pathways are activated vs inhibited - and explain why this is not appropriate for metabolic pathways."

### Prokaryotic / non-model
> "This is a Pseudomonas aeruginosa RNA-seq DE list with PA-locus-tag gene IDs. Run KEGG enrichment with the right organism code and keyType, without forcing an OrgDb or bitr."

### Reproducibility
> "Pin the current human KEGG release as a snapshot, record the date, and run my enrichment against the snapshot so a rerun next year gives the same pathways."

### Multi-condition and modules
> "Compare KEGG enrichment between my up- and down-regulated gene sets in one faceted dotplot, and also run KEGG module enrichment to localize which sub-process is hit."

## What the Agent Will Do
1. Identify the question (membership, ranking, or signed perturbation) and pick enrichKEGG, gseKEGG, or SPIA accordingly.
2. Convert gene IDs to the type KEGG expects (Entrez for eukaryotes, locus tags for prokaryotes) and build the measured universe.
3. Verify the organism code with search_kegg_organism when the organism is not a common model.
4. Run the chosen method with documented thresholds and an explicit universe.
5. For reproducibility, snapshot the KEGG release with gson and record the access date.
6. Translate result IDs to symbols with setReadable (eukaryotes only) and report p.adjust/qvalue with fold enrichment, not raw p-values.
7. Hand plotting to enrichment-visualization, or overlay data on the KEGG map with pathview.

## Common Organism Codes

| Code | Organism | Notes |
|------|----------|-------|
| hsa | Human | Entrez == KEGG gene ID |
| mmu | Mouse | Entrez == KEGG gene ID |
| rno | Rat | Entrez == KEGG gene ID |
| dre | Zebrafish | |
| dme | Drosophila | |
| cel | C. elegans | |
| sce | S. cerevisiae | |
| ath | Arabidopsis | |
| eco | E. coli K-12 | Bacterial; locus tags (b-numbers) |
| pae | P. aeruginosa PAO1 | Bacterial; locus tags (PA-numbers) |
| bsu | B. subtilis 168 | Bacterial |
| ko | KEGG Orthology | Cross-species; use with KO IDs for non-model organisms |

Use `search_kegg_organism('species_name', by = 'scientific_name')` to find codes for other organisms. KEGG covers thousands of species.

## Understanding Results

| Column | Description |
|--------|-------------|
| ID | KEGG pathway/module ID (hsa04110, M00001) |
| Description | Pathway/module name |
| GeneRatio | Query genes in the set / query genes mapped to any set |
| BgRatio | Set genes in universe / universe genes mapped |
| pvalue | Raw p-value |
| p.adjust | BH-adjusted p-value (report this) |
| qvalue | q-value |
| geneID | Genes in the set (raw IDs until setReadable) |
| Count | Number of query genes in the set |

SPIA output adds NDE, pNDE (over-representation), tA and pPERT (perturbation through the topology), pG/pGFdr (combined global), and Status (Activated/Inhibited).

## Tips
- Eukaryotes: KEGG needs Entrez gene IDs. Convert from SYMBOL/ENSEMBL with bitr() and set keyType='ncbi-geneid'.
- Bacteria/prokaryotes: pass locus tags directly with keyType='kegg'; no bitr() and no OrgDb. setReadable() cannot run without an OrgDb.
- Non-model organisms with no KEGG genome: map proteins to KO (eggNOG-mapper, BlastKOALA, KofamScan) and enrich with organism='ko'.
- Always specify the universe (the measured genes); the default background is all KEGG-annotated genes and inflates significance.
- Pin the KEGG release with gson_KEGG() and run enricher()/GSEA() against the snapshot. use_internal_data=TRUE does NOT pin current KEGG; it loads the deprecated 2012 KEGG.db.
- SPIA needs named log2 fold-changes plus the universe and is signaling-only; do not run it on metabolic maps such as glycolysis.
- When comparing conditions, compare pathway-ID sets or use compareCluster(); never compare raw p-values across conditions or across the live vs SPIA-bundled snapshot.
- enrichMKEGG() tests modules (M-numbers): higher resolution, lower power, sparser coverage than full maps.
- Report fold enrichment (GeneRatio / BgRatio) and p.adjust, not raw p-values.
- See enrichment-visualization for dot/cnet/emap plotting; use pathview only for overlaying data on the KEGG map image.

## Related Skills

- go-enrichment - Hypergeometric ORA and the background-universe problem
- gsea - GSEA running-sum engine and ranking-metric choice
- reactome-pathways - Reactome curated-pathway enrichment (reproducible local DB)
- wikipathways - WikiPathways community-pathway enrichment
- enrichment-visualization - Dot/bar/cnet/emap plots of enrichment results
- differential-expression/de-results - Source of the gene list and the fold-changes
- workflows/expression-to-pathways - End-to-end DE-to-enrichment pipeline
