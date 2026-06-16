---
name: bio-pathway-go-enrichment
description: Runs Gene Ontology over-representation analysis (ORA) on a gene LIST with clusterProfiler enrichGO, the one-sided hypergeometric/Fisher 2x2 test phyper(k-1, M, N-M, n, lower.tail=FALSE). Covers why the BACKGROUND universe (not the gene list) is the null and decides significance, why omitting universe= is a bug, why enrichGO defaults to ont='MF' not 'BP', why pvalueCutoff filters p.adjust not raw p, why ORA discards effect magnitude and inherits GO-DAG true-path redundancy (simplify, topGO), why RNA-seq gene-length bias inflates long-gene terms (GOseq Wallenius), plus GeneRatio/BgRatio, bitr ID mapping, minGSSize/maxGSSize, groupGO. Use when a pre-selected gene list (DE hits, co-expression module, screen, GWAS-mapped) needs GO annotation. For a ranked no-cutoff analysis see gsea; for other databases see kegg-pathways, reactome-pathways, wikipathways; DE source is differential-expression/de-results; plots in enrichment-visualization.
tool_type: r
primary_tool: clusterProfiler
---

## Version Compatibility

Reference examples tested with: clusterProfiler 4.18.4+, org.Hs.eg.db 3.22+ (goseq 1.54+ for the length-bias snippet).

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

GO annotation lives in the local org.*.eg.db OrgDb and GO.db, both pinned to the Bioconductor release, so a GO ORA is reproducible given the package versions - record `packageVersion('org.Hs.eg.db')` and `packageVersion('GO.db')` with results. enrichGO moved no core arguments recently, but several plot helpers migrated to enrichplot in clusterProfiler 4.x (those live in enrichment-visualization).

# GO Over-Representation Analysis

**"Which biological processes are enriched in my gene list?"** -> Test each GO term for over-representation of the query genes against a defined background with the one-sided hypergeometric test - because the BACKGROUND universe, not the gene list, is what decides which terms look significant.
- R: `enrichGO(gene, universe, OrgDb, keyType='ENTREZID', ont='BP')`

Scope: hypergeometric ORA of a gene LIST against GO terms, with background-universe selection, ID conversion, GO-DAG redundancy reduction, RNA-seq length-bias correction, and the generic `enricher` test for custom gene sets. A ranked-list / no-cutoff analysis -> gsea. KEGG/Reactome/WikiPathways gene sets -> kegg-pathways, reactome-pathways, wikipathways. The DE list source -> differential-expression/de-results. Plots -> enrichment-visualization.

## The Single Most Important Modern Insight -- ORA Is a Competitive 2x2 Hypergeometric Test Whose Null IS the Chosen Universe

ORA does not answer "which pathways are in my gene list". It is a competitive gene-sampling test (Goeman & Buhlmann 2007 *Bioinformatics* 23:980): of the genes flagged (the foreground), are more annotated to term T than expected when drawing the same number at random from the universe? The p-value is the upper tail of the hypergeometric, computed verbatim by DOSE/clusterProfiler as `phyper(k-1, M, N-M, n, lower.tail=FALSE)` = P(X>=k), the one-sided Fisher exact test on the 2x2 table. Here N = universe genes carrying any GO annotation, M = universe genes in T, n = foreground genes annotated, k = the overlap (the `Count` column). The report columns are GeneRatio = k/n and BgRatio = M/N - both denominators restricted to ANNOTATED genes - and fold enrichment = GeneRatio/BgRatio.

Three consequences drive every misuse:

1. **The universe is the null, not a setting.** Change N or M and every p-value changes. Omitting `universe=` defaults N to ALL annotated genes (~18k for human BP); if the assay only measured ~12k genes, terms for tissue-restricted and lowly-expressed genes go spuriously significant. Omitting `universe=` is a bug, not a default - set it to the genes that COULD have entered the foreground (the tested-gene set), map foreground and universe identically, and report N. The whole-genome background is defensible only when every gene truly could have been detected (Wijesooriya 2022; Timmons 2015).
2. **ORA throws away magnitude and inherits the GO DAG.** A gene is in or out at one threshold (no effect size), so a 2000-gene term at 1.2x fold can beat a 12-gene term at 4x on p-value alone - always read fold enrichment alongside p.adjust. The true-path rule propagates each annotation to all ancestors, so one real signal lights up a whole lineage ("cell cycle", "cell cycle process", "mitotic cell cycle" together); those tests are positively correlated, BH still valid but the term list over-reports. Resolve with `simplify()` (semantic collapse, per ontology) or topGO elim/weight (decorrelation in the test).
3. **The deliverable is never "the enriched pathways."** It is a correctly-backgrounded, effect-sized, redundancy-resolved short list of HYPOTHESES - and it cannot be used to validate the DE list that produced it (circular: the terms are a deterministic function of the same genes; Timmons 2015).

## ORA vs GSEA (the central fork)

ORA needs a pre-selected LIST plus a BACKGROUND and binarizes significant/not; GSEA needs a RANKED vector of ALL genes and no cutoff. Pick by whether a ranking exists and whether the cutoff would be arbitrary. The full three-generations taxonomy (ORA vs FCS vs topology) and competitive-vs-self-contained null theory live in the category README - this skill owns the ORA/GO slice.

| Scenario | Method | Why |
|----------|--------|-----|
| All genes carry a DE statistic, cutoff would be arbitrary | GSEA (gseGO) -> gsea | uses the full ranking; no threshold |
| Pre-selected list (co-expression module, GWAS-mapped, screen hits, markers) | ORA (enrichGO) | no ranking available; ORA is appropriate |
| Very small list (< ~15-20 genes) | low ORA power; report fold enrichment + counts, consider GSEA | hypergeometric power collapses on tiny lists |
| RNA-seq DE list with length/selection bias | GOseq (Wallenius) | length-corrected ORA; standard ORA inflates long-gene terms |

## Tool Taxonomy

| Source / method | Citation | Mechanism / role | When |
|-----------------|----------|------------------|------|
| enrichGO (clusterProfiler) | Yu 2012 *OMICS* 16:284; Wu 2021 *Innovation* 2:100141 | one-sided hypergeometric per GO term; local OrgDb | the default ORA workhorse for a gene list |
| GO DAG (BP/MF/CC) | Ashburner 2000 *Nat Genet* 25:25 | three DAGs; true-path propagation to ancestors | the annotation structure being tested |
| simplify (GOSemSim) | Wang 2007 *Bioinformatics* 23:1274 | semantic-similarity de-redundancy, per ontology | collapse redundant ancestor lineages, keep calibrated p/FDR |
| topGO elim/weight/weight01 | Alexa 2006 *Bioinformatics* 22:1600 | decorrelates the GO graph inside the test | specificity-resolved short list (treat scores as ranking, not FDR) |
| GOseq | Young 2010 *Genome Biol* 11:R14 | Wallenius noncentral hypergeometric weighted by a length PWF | RNA-seq DE with gene-length/selection bias |
| enricher (clusterProfiler) | Yu 2012 *OMICS* 16:284 | same hypergeometric engine on a custom TERM2GENE | any gene set (MSigDB, in-house) not in a DB function |
| gseGO / GSEA | (route -> gsea) | rank-based running-sum, permutation null | a ranking exists; no arbitrary cutoff |

## Run the GO ORA

**Goal:** Find GO terms over-represented in a gene list relative to the genes that could have been selected.

**Approach:** Build the foreground and the universe with the SAME ID mapping, set `ont` explicitly (the source default is 'MF'), pass `universe=` (omitting it is a bug), and read fold enrichment alongside p.adjust.

```r
library(clusterProfiler)
library(org.Hs.eg.db)

ego <- enrichGO(gene          = gene_list,        # foreground ENTREZ IDs
                universe      = universe_ids,     # tested-gene set, mapped identically -- NOT the genome
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = 'BP',             # SET explicitly: source default is 'MF', not 'BP'
                pAdjustMethod = 'BH',
                pvalueCutoff  = 0.05,             # filters p.adjust (despite the name), not raw pvalue
                qvalueCutoff  = 0.2,
                minGSSize     = 10,
                maxGSSize     = 500,
                readable      = TRUE)             # map ENTREZ -> SYMBOL in the output
```

The returned `enrichResult` has columns `ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count` (plus `ONTOLOGY` when `ont='ALL'`). `pvalueCutoff` filters the ADJUSTED p, so an empty table usually means the cutoff or the universe, not biology - inspect everything with `pvalueCutoff=1, qvalueCutoff=1`.

## Build the Foreground and Universe from DE Results

**Goal:** Turn a DE table into the foreground gene vector and the matched background universe.

**Approach:** Filter the DE table to the hits for the foreground; take the genes that were actually TESTED for the universe (DESeq2: rows with non-NA pvalue survive independent filtering); map both with the same `bitr` call. The DE mechanics and the `$padj`/`$adj.P.Val` column choice live in differential-expression/de-results.

```r
de <- read.csv('de_results.csv')

sig_genes  <- de$gene_id[de$padj < 0.05 & abs(de$log2FoldChange) > 1]   # foreground = hits
all_tested <- de$gene_id[!is.na(de$pvalue)]                            # universe = tested genes, NOT all rows, NOT the genome

fg_map <- bitr(sig_genes,  fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
bg_map <- bitr(all_tested, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

gene_list    <- unique(fg_map$ENTREZID)   # deduplicate one-to-many maps before counting
universe_ids <- unique(bg_map$ENTREZID)
```

`bitr` one-to-many maps produce duplicate rows that inflate `Count`; deduplicate. If more than ~15% of genes fail to convert the result is unreliable - report the conversion rate. Mixed up- and down-regulated genes cancel in one list: run ORA separately per direction when direction matters.

## Reduce GO-DAG Redundancy with simplify

**Goal:** Collapse the redundant ancestor lineage so one biological signal is one entry, not a dozen.

**Approach:** `simplify()` removes terms whose semantic similarity to a kept term exceeds the cutoff. It operates on ONE ontology (GOSemSim defines similarity within a single DAG), so run BP/MF/CC separately and simplify each - it does NOT de-redundify an `ont='ALL'` object.

```r
ego_bp <- enrichGO(gene_list, universe = universe_ids, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = 'BP', readable = TRUE)
ego_bp <- simplify(ego_bp, cutoff = 0.7, by = 'p.adjust', select_fun = min, measure = 'Wang')
```

`measure='Wang'` (the default) is graph-topology-based and stable across annotation releases; IC-based measures ('Resnik', 'Lin', 'Jiang', 'Rel') shift with the annotation corpus. topGO elim/weight01 is the alternative that decorrelates inside the test, returning a specificity-resolved list directly - but its conditioned p-values are best treated as a ranking, not calibrated FDR (Alexa 2006).

## Correct RNA-seq Length Bias with GOseq

**Goal:** Stop long, highly-expressed genes from looking enriched for a purely technical reason.

**Approach:** DE-detection power scales with read count, which scales with transcript length and expression, so the foreground is enriched for long genes - and RPKM/TMM normalization does NOT fix it (it corrects abundance, not detection power). GOseq fits a probability weighting function (PWF) over the bias variable and tests with the Wallenius noncentral hypergeometric (Young 2010). The input is a NAMED 0/1 vector over ALL tested genes; goseq returns UNADJUSTED p-values, so apply BH afterward.

```r
library(goseq)

all_genes <- de$gene_id[!is.na(de$pvalue)]
de_genes  <- as.integer(all_genes %in% sig_genes)   # named binary vector over the tested set
names(de_genes) <- all_genes

pwf <- nullp(de_genes, 'hg38', 'ensGene')           # fits the length PWF; inspect the fit plot
go  <- goseq(pwf, 'hg38', 'ensGene', method = 'Wallenius')   # default; 'Hypergeometric' ignores bias (= standard ORA)
go$padj <- p.adjust(go$over_represented_pvalue, method = 'BH')   # goseq does NOT BH-correct internally
```

GSEA on a length-neutral ranking statistic (the moderated t / Wald z) is largely immune to this bias - one more reason to consider gsea for RNA-seq.

## All Three Ontologies and a Descriptive Breakdown

`ont='ALL'` runs BP/MF/CC separately and rbinds them with an `ONTOLOGY` column (`pool=FALSE` default; `pool=TRUE` treats the three as one set). `groupGO` is NOT a test - it classifies genes at a fixed DAG level for a GO-slim overview (counts, no p-values); never read its counts as significance.

```r
ego_all <- enrichGO(gene_list, universe = universe_ids, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = 'ALL', readable = TRUE)
ggo     <- groupGO(gene_list, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = 'BP', level = 3, readable = TRUE)
```

## Custom Gene Sets with enricher

For gene sets not covered by a DB function (MSigDB collections, in-house sets), `enricher` runs the SAME hypergeometric engine against a two-column TERM2GENE table; pass the same explicit `universe`.

```r
ego_custom <- enricher(gene_list, TERM2GENE = t2g, universe = universe_ids,
                       pvalueCutoff = 0.05, pAdjustMethod = 'BH', minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2)
```

## Other Organisms

Swap the OrgDb: `org.Mm.eg.db` (mouse), `org.Dr.eg.db` (zebrafish), `org.Sc.sgd.db` (yeast, `keyType='ORF'`). Check usable key types with `keytypes(OrgDb)`.

## Per-Method Failure Modes

### Whole-genome or default universe
**Trigger:** omitting `universe=`, or passing the genome when the assay measured fewer genes. **Mechanism:** N defaults to all annotated genes, inflating the denominator with genes that never could have been selected. **Symptom:** a confident table where tissue-restricted / lowly-expressed-gene terms dominate. **Fix:** set `universe=` to the tested-gene set, map foreground and universe identically, report N.

### p read without fold enrichment (term-size trap)
**Trigger:** ranking results by p.adjust alone. **Mechanism:** a 2000-gene term has enormous power at tiny fold enrichment; p scales with term size. **Symptom:** vague broad terms ("cellular process") top the list, specific terms buried. **Fix:** read fold enrichment = (k/n)/(M/N) alongside p.adjust; trim extremes with minGSSize=10, maxGSSize=500.

### Redundant ancestor lineage counted as findings
**Trigger:** reporting "cell cycle", "cell cycle process", "mitotic cell cycle" as separate discoveries. **Mechanism:** true-path propagation lights up a whole lineage from one signal; the tests are positively correlated. **Symptom:** the top 20 is one biological theme repeated. **Fix:** `simplify()` per ontology, or topGO weight01; never count lineage members as independent hits.

### RNA-seq length / selection bias
**Trigger:** standard ORA on an RNA-seq DE list without length correction. **Mechanism:** detection power scales with count ~ length/expression; TMM/RPKM fixes abundance, not power. **Symptom:** long-gene categories (ECM, adhesion) enriched, short-gene (ribosomal) depleted - and it survives FDR. **Fix:** GOseq with a length PWF + `method='Wallenius'`, then BH; or GSEA on a bias-neutral statistic.

### Wrong ID type or silent gene loss
**Trigger:** passing ENSEMBL/SYMBOL with a mismatched `keyType`, or not checking the bitr conversion rate. **Mechanism:** unmapped IDs are dropped, shrinking the foreground; one-to-many maps inflate Count. **Symptom:** "no gene can be mapped", or a suspiciously small/large Count. **Fix:** match `keyType` to one of `keytypes(OrgDb)`, deduplicate after bitr, report conversion rate (flag >15% loss).

### pvalueCutoff misread as raw-p filter
**Trigger:** concluding "no significant terms" when strong raw p exists. **Mechanism:** `pvalueCutoff` filters p.adjust, not pvalue. **Symptom:** an empty table despite plausible signal. **Fix:** inspect with `pvalueCutoff=1, qvalueCutoff=1`, then judge on p.adjust.

### simplify on ont='ALL'
**Trigger:** calling `simplify()` on an `ont='ALL'` object. **Mechanism:** semantic similarity is defined within ONE ontology, not across BP/MF/CC. **Symptom:** redundancy not removed, or an error. **Fix:** run BP/MF/CC separately and simplify each.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `pvalueCutoff = 0.05` | clusterProfiler default | filters on p.adjust (NOT raw pvalue); standard FDR gate |
| `qvalueCutoff = 0.2` | clusterProfiler default | secondary q-value gate; loosen to 1 to inspect all terms |
| `pAdjustMethod = 'BH'` | Benjamini-Hochberg | controls FDR; valid under the positive dependence of true-path-correlated terms (Bonferroni is needlessly strict here) |
| `minGSSize = 10` | enrichGO default | drop tiny sets that overfit and are noisy |
| `maxGSSize = 500` | enrichGO default | drop huge general sets that always "enrich" with trivial fold |
| `simplify(cutoff = 0.7)` | GOSemSim/Wang | semantic-similarity redundancy cutoff; lower keeps more terms, higher is more aggressive |
| fold enrichment > 2 | heuristic | (k/n)/(M/N); a rough "strong" flag, never a substitute for p.adjust |
| ID-conversion loss > 15% | heuristic | above this the foreground is too eroded to trust; report the rate |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `--> No gene can be mapped` | wrong keyType / OrgDb, or IDs not in the OrgDb | match keyType to `keytypes(OrgDb)`; bitr to ENTREZID first |
| Empty result table | `pvalueCutoff` filters p.adjust; or universe too large; or IDs lost | set cutoffs to 1 to inspect; fix the universe; check conversion rate |
| Vague broad terms dominate | ranking by p alone (term-size trap) | read fold enrichment; trim with minGSSize/maxGSSize |
| Many redundant ancestor terms | GO-DAG true-path propagation | `simplify()` per ontology, or topGO weight01 |
| simplify does nothing / errors on ALL | similarity is per-ontology | run BP/MF/CC separately |
| Description column shows IDs not names | not readable | `readable=TRUE` or `setReadable(ego, OrgDb, 'ENTREZID')` |
| Tested MF when expecting BP | enrichGO default `ont='MF'` | set `ont` explicitly every call |

## References

- Ashburner M, Ball CA, Blake JA, et al. 2000. Gene Ontology: tool for the unification of biology. *Nat Genet* 25:25-29.
- Yu G, Wang LG, Han Y, He QY. 2012. clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS* 16:284-287.
- Wu T, Hu E, Xu S, et al. 2021. clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. *The Innovation* 2(3):100141.
- Goeman JJ, Buhlmann P. 2007. Analyzing gene expression data in terms of gene sets: methodological issues. *Bioinformatics* 23:980-987.
- Alexa A, Rahnenfuhrer J, Lengauer T. 2006. Improved scoring of functional groups from gene expression data by decorrelating GO graph structure. *Bioinformatics* 22(13):1600-1607.
- Young MD, Wakefield MJ, Smyth GK, Oshlack A. 2010. Gene ontology analysis for RNA-seq: accounting for selection bias. *Genome Biol* 11(2):R14.
- Wang JZ, Du Z, Payattakool R, et al. 2007. A new method to measure the semantic similarity of GO terms. *Bioinformatics* 23(10):1274-1281.
- Timmons JA, Szkop KJ, Gallagher IJ. 2015. Multiple sources of bias confound functional enrichment analysis of global -omics data. *Genome Biol* 16:186.
- Wijesooriya K, Jadaan SA, Perera KL, et al. 2022. Urgent need for consistent standards in functional enrichment analysis. *PLoS Comput Biol* 18(3):e1009935.

## Related Skills

- gsea - Ranked-list GSEA alternative when a full ranking exists and a cutoff is arbitrary
- kegg-pathways - KEGG pathway and module enrichment
- reactome-pathways - Reactome curated-pathway enrichment
- wikipathways - WikiPathways community-pathway enrichment
- enrichment-visualization - Dot/bar/cnet/emap/tree plots of enrichment results
- differential-expression/de-results - Source of the gene list and the tested-gene universe
- database-access/entrez-fetch - Fetch gene annotations / ID maps from NCBI
- workflows/expression-to-pathways - End-to-end DE-to-enrichment pipeline
