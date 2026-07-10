---
name: bio-workflows-expression-to-pathways
description: 'Orchestrates the full path from differential expression results to redundancy-collapsed functional enrichment: choose ORA vs GSEA, convert gene IDs per method, run enrichGO/enrichKEGG/enrichPathway/enrichWP or gseGO/gseKEGG (clusterProfiler, ReactomePA, rWikiPathways), and visualize. Use when a DESeq2/edgeR/limma result must become enriched GO terms, KEGG/Reactome/WikiPathways pathways, or a GSEA leading edge; when the input is a full ranking for all genes (GSEA, named decreasing vector) or only a pre-selected list (ORA plus a defensible background universe); or when assembling DE-to-pathway end to end. The DE list and ranking statistic come from differential-expression/de-results; per-method nuance lives in the pathway-analysis skills.'
tool_type: r
primary_tool: clusterProfiler
workflow: true
depends_on:
  - pathway-analysis/go-enrichment
  - pathway-analysis/gsea
  - pathway-analysis/kegg-pathways
  - pathway-analysis/reactome-pathways
  - pathway-analysis/wikipathways
  - pathway-analysis/enrichment-visualization
qc_checkpoints:
  - input_validation: "Gene IDs match the method (OrgDb keyType / kegg-id / ENTREZ); >85% convert; background = testable genes"
  - generation_choice: "ORA-vs-GSEA fork decided BEFORE running; a ranking for all genes -> GSEA, a pre-selected list -> ORA"
  - reproducibility: "Tool + database version/date, ranking metric, p-adjust method, and universe recorded; set.seed for GSEA"
---

## Version Compatibility

Reference examples tested with: clusterProfiler 4.10+, org.Hs.eg.db 3.18+, ReactomePA 1.46+, enrichplot 1.22+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Expression to Pathways Workflow

**"Find enriched pathways from my differential expression results"** -> Decide the generation (ORA vs GSEA) FIRST, then convert IDs to the form each method needs, run enrichment against the chosen database, and collapse redundancy before interpreting - because the enrichment result is a claim conditioned on the method, the background universe, and the database version, not a discovery the algorithm hands back.
- R: `enrichGO(...)` / `gseGO(...)` / `enrichKEGG(...)` / `enrichPathway(...)` (clusterProfiler, ReactomePA)

Scope: the ORCHESTRATION of a DE-to-enrichment pipeline - the generation fork, per-method ID conversion, the universe decision, the live-vs-local database caveat, and the handoff to redundancy-collapsed visualization. This workflow does NOT re-teach each method. The null/universe/reproducibility theory and the master method-selection tree -> the pathway-analysis README and the per-method skills (go-enrichment, gsea). ORA mechanics -> go-enrichment; GSEA mechanics -> gsea; per-database IDs and live-DB behavior -> kegg-pathways, reactome-pathways, wikipathways; the DE list and ranking statistic -> differential-expression/de-results; plotting -> enrichment-visualization.

## The Single Most Important Modern Insight -- The First Decision Is the Generation, and It Is Set by the Input, Not by Preference

Pathway analysis has three generations (Khatri 2012 *PLoS Comput Biol* 8:e1002375): over-representation analysis (ORA), functional class scoring / GSEA (FCS), and pathway topology. A workflow that "runs enrichment" without first deciding which generation applies has already made the choice silently - usually ORA, the worst-calibrated corner of the space. The fork is mechanical:

1. **Is there a meaningful per-gene ranking for (nearly) ALL measured genes?** A signed test statistic, the DESeq2 Wald `stat`, or `-sign(log2FC)*log10(p)` for every tested gene -> **GSEA** (a NAMED vector sorted in DECREASING order; the ranking metric IS the experiment). No arbitrary cutoff; detects coordinated weak shifts that ORA misses.
2. **Only a pre-selected LIST (DE hits past a cutoff, a co-expression module, GWAS loci, screen hits)?** -> **ORA**, and the deliverable hinges on a defensible **background universe** - the genes that were testable, not the whole genome. ORA's p-value is whatever the denominator says it is.

The dangerous default is running ORA on data that has a full ranking (binarizing away the signal) or running ORA against the genome (measuring expression, not enrichment). Decide the fork out loud, record it, and record the why (see the pathway-analysis README) - this workflow owns the routing, not the derivation.

## Pipeline Overview

```
DE results (differential-expression/de-results)
    |
    v
[0. Decide the generation: ranking for all genes? -> GSEA | pre-selected list? -> ORA]
    |
    +--> ORA branch: define the TESTABLE-gene universe, convert IDs per method
    |        +--> enrichGO     (OrgDb keyType)        -> go-enrichment
    |        +--> enrichKEGG   ('kegg' / 'ncbi-geneid', LIVE DB)  -> kegg-pathways
    |        +--> enrichPathway (ENTREZ, local DB)    -> reactome-pathways
    |        +--> enrichWP      (ENTREZ, LIVE GMT)    -> wikipathways
    |
    +--> GSEA branch: build a NAMED decreasing vector of ALL genes, set.seed
    |        +--> gseGO / gseKEGG / GSEA(+msigdbr)    -> gsea
    |
    v
[Redundancy collapse + visualization: simplify, pairwise_termsim, dotplot/emapplot/gseaplot2]   (enrichment-visualization)
    |
    v
A claim conditioned on universe + method + database version (record provenance)
```

## Stage Map

| Stage | Goal | Owns the nuance |
|-------|------|-----------------|
| 0. Decide generation | ORA vs GSEA from the available input | pathway-analysis README (method selection) |
| 1. Prepare input | Build the gene list AND/OR the named ranked vector; define the universe | differential-expression/de-results (the stat); pathway-analysis/go-enrichment (the universe rule) |
| 2. Convert IDs | Map to the form each method needs (OrgDb keyType / kegg-id / ENTREZ) | go-enrichment, kegg-pathways, reactome-pathways |
| 3a. ORA | Hypergeometric test of the list vs background | go-enrichment, kegg-pathways, reactome-pathways, wikipathways |
| 3b. GSEA | Running-sum over the full ranking | gsea |
| 4. Collapse + visualize | Reduce redundancy, then plot | enrichment-visualization |

## Decision Tree by Scenario

| Scenario | Route | Why |
|----------|-------|-----|
| All genes carry a DE statistic, a cutoff would be arbitrary | GSEA (gseGO/gseKEGG) -> gsea | uses the full ranking; no cutoff; named decreasing vector |
| Pre-selected list (module, GWAS loci, screen hits), no ranking | ORA (enrichGO/enrichKEGG) -> go-enrichment | no ranking available; define the universe |
| Broad function annotation | enrichGO / gseGO -> go-enrichment, gsea | GO is the broadest LOCAL resource (reproducible) |
| Metabolic / signaling pathways | enrichKEGG / gseKEGG -> kegg-pathways | KEGG maps query a LIVE DB (pin the date) |
| Reaction-level, peer-reviewed, reproducible offline | enrichPathway -> reactome-pathways | local reactome.db, version-pinned |
| Disease/drug sets the others miss, broad species | enrichWP -> wikipathways | community-curated; LIVE versioned GMT |
| Bacterial / prokaryotic data | enrichKEGG with locus tags + KEGG organism code -> kegg-pathways | KEGG covers prokaryotes; OrgDb usually does not |
| RNA-seq with strong gene-length bias | GOseq -> go-enrichment | length-aware ORA null |
| Multiple conditions/clusters side by side | compareCluster -> any DB | one model, faceted dotplot; never compare p across separate runs |
| The DE list / ranking statistic itself | -> differential-expression/de-results | that is upstream, not enrichment |
| Why this null, which universe, version reporting | -> go-enrichment (universe), gsea (null) | per-method theory owned by each skill |

## Stage 1: Prepare the Input (list, ranked vector, and the universe)

**Goal:** Turn a DE table into the two possible inputs - a gene LIST for ORA and a NAMED decreasing vector for GSEA - and define the background universe as the testable genes.

**Approach:** Read the DE result, derive the significant list, build the ranked vector from a signed statistic (not a bare log2FC), and set the universe to exactly the genes that entered the DE test. The DE mechanics (the `$padj` vs `$adj.P.Val` column, shrinkage) live at differential-expression/de-results - this is only input shaping.

```r
library(clusterProfiler)
library(org.Hs.eg.db)

res <- read.csv('deseq2_results.csv', row.names = 1)

# ORA input: a pre-selected list (DESeq2 padj column; limma/edgeR name it differently)
sig_genes <- rownames(subset(res, padj < 0.05 & abs(log2FoldChange) > 1))

# Background universe = genes that were TESTABLE (entered the DE test), NOT the genome.
# Using the genome measures expression bias, not enrichment (the universe rule; see pathway-analysis/go-enrichment).
universe_genes <- rownames(res[!is.na(res$pvalue), ])

# GSEA input: a NAMED vector of ALL genes, sorted DECREASING by a signed metric.
# Prefer the Wald stat (magnitude + precision); a bare log2FC over-weights noisy low-count genes.
# TRAP: a table from lfcShrink(type='apeglm'/'ashr') has NO `stat` column -- shrinkage drops it.
# Rank from the UNSHRUNK results(dds)$stat; use edgeR `sign(logFC)*-log10(PValue)` or limma `t`.
ranked <- res$stat
names(ranked) <- rownames(res)
ranked <- sort(ranked[!is.na(ranked)], decreasing = TRUE)
```

## Stage 2: Convert Gene IDs to What Each Method Needs

**Goal:** Map identifiers to the exact ID type each enrichment function expects, because a mismatch returns zero hits silently.

**Approach:** Use `bitr` (OrgDb) for SYMBOL/ENSEMBL -> ENTREZ, keep both list and ranked vector in the same ID space, deduplicate, and track the conversion rate. Per-method ID rules are owned by each DB skill; the table below is the routing summary.

```r
# enrichGO accepts ENSEMBL/SYMBOL/ENTREZ via keyType=; ENTREZ is the safe lingua franca downstream
sig_entrez <- bitr(sig_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
bg_entrez <- bitr(universe_genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

# Carry the ranking through conversion: name the kept stat by its ENTREZ id
ranked_map <- bitr(names(ranked), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
ranked_list <- ranked[ranked_map$SYMBOL]
names(ranked_list) <- ranked_map$ENTREZID
ranked_list <- ranked_list[!duplicated(names(ranked_list))]   # dedup or GSEA biases the score
ranked_list <- sort(ranked_list, decreasing = TRUE)           # re-sort: bitr remap can reorder rows

conv_rate <- nrow(sig_entrez) / length(sig_genes)   # report it; <0.85 -> wrong ID type/organism
```

| Method | keyType / ID required | Convert with |
|--------|-----------------------|--------------|
| enrichGO / gseGO | OrgDb keyType ('ENSEMBL', 'SYMBOL', 'ENTREZID') | bitr |
| enrichKEGG / gseKEGG | 'kegg' or 'ncbi-geneid' (NOT ENSEMBL/OrgDb) | bitr to ENTREZID, pass keyType='ncbi-geneid' (bitr_kegg only converts among KEGG ID flavors) |
| enrichPathway / gsePathway (ReactomePA) | ENTREZ | bitr |
| enrichWP / gseWP (WikiPathways) | ENTREZ + organism string | bitr |

## Stage 3a: ORA Branch (pre-selected list + universe)

**Goal:** Test each gene set for over-representation of the list against the testable-gene background.

**Approach:** Always pass `universe=`; run GO ontologies separately; KEGG/Reactome/WikiPathways each need their own ID form. KEGG and WikiPathways query a LIVE database (internet-dependent, not reproducible across releases - pin the run date); GO and Reactome read local annotation (reproducible given the Bioconductor release).

```r
# GO ORA - universe is the decision; simplify() collapses DAG redundancy (BP/MF/CC separately, not 'ALL')
go_bp <- enrichGO(sig_entrez$ENTREZID, universe = bg_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                  ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, readable = TRUE)
go_bp <- simplify(go_bp, cutoff = 0.7, by = 'p.adjust')

# KEGG ORA - LIVE KEGG REST API; needs internet; record the access date for reproducibility
kegg <- enrichKEGG(sig_entrez$ENTREZID, universe = bg_entrez$ENTREZID, organism = 'hsa', keyType = 'ncbi-geneid', pvalueCutoff = 0.05)   # Entrez input; keyType='kegg' = Entrez for eukaryotes / locus tags for prokaryotes, 'ncbi-geneid' is explicit
kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')

# Reactome ORA - ENTREZ required; LOCAL reactome.db so reproducible given the release
library(ReactomePA)
reactome <- enrichPathway(sig_entrez$ENTREZID, universe = bg_entrez$ENTREZID, organism = 'human', pvalueCutoff = 0.05, readable = TRUE)
```

## Stage 3b: GSEA Branch (named decreasing vector of all genes)

**Goal:** Find gene sets whose genes shift coordinately across the full ranking, without a significance cutoff.

**Approach:** Run on the named decreasing `ranked_list`, fix the permutation seed so p-values are reproducible, then read the leading edge as the interpretable core. clusterProfiler GSEA is preranked / gene-permutation (the inter-gene-correlation-UNcorrected null) - a discovery screen; see pathway-analysis/gsea for the calibration caveat (CAMERA/ROAST).

```r
set.seed(123)   # permutation reproducibility; without it p-values drift across runs

gsea_go <- gseGO(ranked_list, OrgDb = org.Hs.eg.db, ont = 'BP',
                 minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)

gsea_kegg <- gseKEGG(ranked_list, organism = 'hsa',
                     minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
```

## Stage 4: Collapse Redundancy, Then Visualize

**Goal:** Reduce overlapping terms to distinct findings before drawing conclusions, then plot deliberately.

**Approach:** A list of 40 significant GO terms is often a few biological stories told many times (shared genes via the GO true-path rule). Collapse with `simplify`/`pairwise_termsim`, then plot - `emapplot`/`treeplot` require `pairwise_termsim()` first (cnetplot does NOT; it draws the gene-concept network from the `geneID` column directly), and `gseaplot2` is for a gseaResult not an enrichResult. Encoding choice and required pre-steps are owned by pathway-analysis/enrichment-visualization.

```r
library(enrichplot)

go_bp <- pairwise_termsim(go_bp)              # required before emapplot/treeplot
dotplot(go_bp, showCategory = 20)             # GeneRatio vs Count: pick the encoding deliberately
emapplot(go_bp, showCategory = 30)            # redundancy-collapsed term-similarity map
gseaplot2(gsea_go, geneSetID = 1:3)           # gseaResult only, not enrichResult
```

## Multi-Condition Comparison

**Goal:** Compare enrichment across conditions in one model instead of comparing p-values from separate runs.

**Approach:** `compareCluster` fits all gene lists together and facets the dotplot; never compare raw -log10(p) across separate enrichments (it scales with set size and sample size). For GSEA, compare NES, not p.

```r
gene_clusters <- list(A = sig_A, B = sig_B, C = sig_C)
cc <- compareCluster(gene_clusters, fun = 'enrichKEGG', organism = 'hsa',
                     universe = bg_entrez$ENTREZID)   # compareCluster forwards ... to fun; omitting universe silently reverts to the whole-genome background
dotplot(cc, showCategory = 10)
```

## Per-Method Failure Modes

### Whole-genome background
**Trigger:** `universe=` left at default while only ~12k genes were expressed. **Mechanism:** the hypergeometric p-value is fully determined by the denominator; the genome inflates any set whose members are expressed in the tissue. **Symptom:** many tissue-specific terms enrich with tiny p. **Fix:** set `universe` to the genes that entered the DE test (the testable set).

### ORA on a ranked dataset
**Trigger:** filtering all-gene DE results to a list and running ORA. **Mechanism:** binarizing at an arbitrary cutoff discards magnitude and the coordinated-weak signal. **Symptom:** GSEA finds sets ORA missed. **Fix:** if a ranking exists for all genes, run GSEA; reserve ORA for genuinely unranked lists.

### Wrong ID type for the database
**Trigger:** ENSEMBL/SYMBOL passed to enrichKEGG/enrichPathway/enrichWP. **Mechanism:** those expect kegg-id/ENTREZ; unmatched IDs are dropped. **Symptom:** zero terms, no error. **Fix:** `bitr`/`bitr_kegg` to the required ID; check `conv_rate`.

### GSEA without set.seed or with an unsorted vector
**Trigger:** no seed, or a list that is not named and decreasing. **Mechanism:** permutation p-values drift run to run; an unsorted/unnamed vector errors or mis-ranks. **Symptom:** different leading edges each run, or a names error. **Fix:** build the named decreasing vector and `set.seed`.

### Ranking a GSEA vector off a shrunken-LFC object
**Trigger:** building the ranked vector from a `lfcShrink(type='apeglm'/'ashr')` table, or ranking by bare `log2FoldChange`. **Mechanism:** apeglm/ashr DROP the `stat` column, so the vector silently falls back to shrunken LFC; low-count genes with unstable large FC then hijack the leading edge. **Symptom:** the leading edge is dominated by low-baseMean genes, or `res$stat` is NULL. **Fix:** rank from the unshrunken `results(dds)$stat` (DESeq2), limma `topTable$t`, or edgeR `sign(logFC)*-log10(PValue)`; reserve shrunken LFC for visualization.

### Live-DB result reported as reproducible
**Trigger:** KEGG/WikiPathways result with no recorded date. **Mechanism:** those query the current data release; the same code returns different pathways later. **Symptom:** a collaborator cannot reproduce the figure. **Fix:** record the access date and data version; prefer local GO/Reactome when reproducibility is paramount.

### Redundancy read as replication
**Trigger:** interpreting 40 overlapping GO terms as 40 findings. **Mechanism:** the true-path rule and pathway overlap mean shared genes drive many sets. **Symptom:** the same 3-5 genes explain the top 20 terms. **Fix:** `simplify`/`pairwise_termsim`, inspect the leading-edge/`geneID` core, report clusters of terms.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `pvalueCutoff = 0.05` | clusterProfiler default | filters on p.adjust by default in enrichResult; standard FDR gate |
| `qvalueCutoff = 0.2` | clusterProfiler default | secondary q-value gate |
| `pAdjustMethod = 'BH'` | Benjamini-Hochberg | valid FDR control under positive dependence (overlapping sets); Bonferroni over-corrects |
| `minGSSize = 10` | enrichGO/gseGO default | drop tiny sets that overfit |
| `maxGSSize = 500` | enrichGO/gseGO default | drop overly broad sets that always "enrich" |
| `simplify(cutoff = 0.7)` | GOSemSim semantic similarity | GO DAG redundancy cutoff; lower keeps more terms |
| conversion rate > 0.85 | practical QC | <85% ID conversion flags a wrong ID type/organism |
| `set.seed(123)` | reproducibility | any fixed seed; the point is to fix the permutation draw |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| enrichKEGG returns 0 terms | ENSEMBL passed (needs kegg-id/ENTREZ), wrong organism code, or KEGG API down | convert with bitr_kegg; check organism; retry (live DB) |
| `--> No gene can be mapped` | wrong keyType/OrgDb for the input IDs | match keyType to the actual ID type |
| gseGO error about names | vector not named or not sorted decreasing | build a named vector sorted `decreasing = TRUE` |
| emapplot/treeplot empty or errors | `pairwise_termsim()` not run first (cnetplot does not need it) | run `pairwise_termsim()` before emapplot/treeplot |
| simplify fails on ont='ALL' | simplify needs one ontology | run BP/MF/CC separately, then simplify each |
| different results each run | no set.seed, or the live KEGG/WP DB changed | set.seed; pin and record the DB version/date |
| all terms have NA Description | `readable`/`setReadable` not applied | set `readable = TRUE` or call `setReadable` |

## References

- Khatri P, Sirota M, Butte AJ. 2012. Ten years of pathway analysis: current approaches and outstanding challenges. *PLoS Comput Biol* 8:e1002375.
- Subramanian A, Tamayo P, Mootha VK, et al. 2005. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 102:15545-15550.
- Goeman JJ, Buhlmann P. 2007. Analyzing gene expression data in terms of gene sets: methodological issues. *Bioinformatics* 23:980-987.
- Wu T, Hu E, Xu S, et al. 2021. clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. *The Innovation* 2:100141.
- Yu G, He QY. 2016. ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. *Mol BioSyst* 12:477-479.
- Kanehisa M, Goto S. 2000. KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Res* 28:27-30.
- Young MD, Wakefield MJ, Smyth GK, Oshlack A. 2010. Gene ontology analysis for RNA-seq: accounting for selection bias. *Genome Biol* 11:R14.
- Wijesooriya K, Jadaan SA, Perera KL, Kaur T, Ziemann M. 2022. Urgent need for consistent standards in functional enrichment analysis. *PLoS Comput Biol* 18:e1009935.

## Related Skills

- pathway-analysis/go-enrichment - GO over-representation, background universe, redundancy reduction, length bias
- pathway-analysis/gsea - Ranked-list GSEA, named decreasing vector, ranking metric, leading edge, NES
- pathway-analysis/kegg-pathways - KEGG pathway/module enrichment, live DB, prokaryotes, multi-condition
- pathway-analysis/reactome-pathways - Reactome curated-pathway ORA and GSEA, ENTREZ IDs, reproducible local DB
- pathway-analysis/wikipathways - WikiPathways community-pathway enrichment, versioned GMT, broad species
- pathway-analysis/enrichment-visualization - Dot/bar/cnet/emap/GSEA plots and required pre-steps
- differential-expression/de-results - Source of the gene list and the ranking statistic
