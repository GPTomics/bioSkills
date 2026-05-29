---
name: bio-differential-expression-de-results
description: Extracts, filters, annotates, and exports differential expression results from DESeq2 or edgeR with proper handling of padj=NA (independent filtering, Cook's outliers, all-zero), multiple-testing correction choice (BH vs Storey q-value vs IHW vs lfsr), TREAT vs post-hoc fold-change filtering, p-value histogram diagnostics, gene annotation via org.db/biomaRt/mygene, GSEA preranked input, ORA background construction, replication reality (Schurch 2016 small-n result), and SABV/sex-stratified reporting. Use when extracting and interpreting DE results, troubleshooting padj=NA, choosing FDR method, preparing ranked lists for pathway analysis, annotating gene IDs, or comparing DESeq2 vs edgeR outputs.
tool_type: r
primary_tool: DESeq2
---

## Version Compatibility

Reference examples tested with: DESeq2 1.42+, edgeR 4.0+, IHW 1.34+, qvalue 2.34+, ashr 2.2+, AnnotationDbi 1.66+, org.Hs.eg.db 3.18+, biomaRt 2.58+, mygene 1.38+ (Python), dplyr 1.1+, openxlsx 4.2+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# DE Results

**"What are my significant genes?"** -> Extract DE estimates and p-values from the fitted model, handle missing padj correctly, apply FDR control appropriate to the design, and produce the table or ranked list the downstream tool actually needs.

## The Single Most Important Modern Insight -- `padj = NA` has three distinct meanings

A `NA` in the `padj` column is not a missing value; it is a flag indicating which filter excluded the gene. The three causes -- independent filtering, Cook's distance outlier, and all-zero in a group -- have completely different remediations. Dropping all NA rows blindly silently discards real signal, most often from low-count master regulators (transcription factors expressed at ~10 counts) that pass biology but fail the data-driven baseMean threshold.

| `padj = NA` cause | DESeq2 detection | What it means | Fix if undesired |
|-------------------|------------------|---------------|-------------------|
| Independent filtering | finite `pvalue`, `NA` `padj`, baseMean below auto threshold | Removed before BH adjustment to maximize rejections at `alpha` | `results(dds, independentFiltering = FALSE)` OR `filterFun = ihw` |
| Cook's distance outlier | `NA` `pvalue`, `NA` `padj`, baseMean > 0, group has >=3 reps | One sample has Cook's > `qf(0.99, p, m-p)` | `results(dds, cooksCutoff = FALSE)` |
| All-zero or near-zero in a group | `NA` `pvalue` AND baseMean very low | Insufficient information to test | Filter at preprocess time; or accept |

Independent filtering (Bourgon, Gentleman, Huber 2010 *PNAS* 107:9546) chooses the baseMean threshold to maximize rejections. The filter MUST be independent of the test statistic under the null -- this is why baseMean (the across-sample mean) is the canonical choice. Using "min count in treatment group" as a filter VIOLATES the independence requirement and inflates type-I error. Most pipelines unknowingly do this; do not.

A second axis: at n>=7 per group, `DESeq()` also REPLACES outlier counts via `replaceOutliers()` and refits (default `minReplicatesForReplace = 7`). Cook's filtering is NOT computed for continuous covariates -- a continuous-covariate analysis has effectively no outlier filtering.

## Algorithmic Taxonomy

| Method | What it computes | When to use | Failure mode |
|--------|------------------|-------------|--------------|
| BH (`p.adjust(method='BH')`, DESeq2 default `pAdjustMethod='BH'`) | FDR at fixed alpha; Benjamini-Hochberg 1995 | Default for most RNA-seq DE | Assumes independence or PRDS; many overlapping tests violate |
| Storey q-value (`qvalue::qvalue`) | q-value using estimated pi_0 | Genome-scale with many true nulls | Pi_0 estimation can fail at small test counts |
| IHW (`results(filterFun=ihw)`, Ignatiadis 2016 *Nat Methods* 13:577) | Weighted BH with covariate-informed weights | Modern default for DESeq2; +5-20% discoveries at same FDR | Covariate MUST be independent under null |
| ashr local false sign rate (`lfsr` / `svalue`, with `svalue=TRUE`) | P(sign of estimate is wrong) | When effect-direction certainty is what matters | Conservative lower bound on FDR; not interchangeable with padj |
| BY (`p.adjust(method='BY')`) | Benjamini-Yekutieli; arbitrary dependence | Strongly correlated tests | Uniformly conservative; rarely needed for DE |
| Holm / Bonferroni | FWER | Small confirmatory test sets | Far too conservative for genome-scale |
| TREAT / `lfcThreshold=` | FDR for "|LFC| > tau" hypothesis | Pre-specified biologically meaningful threshold | tau must be set BEFORE looking at data |

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Standard bulk DE, two groups | DESeq2 results with default BH; report padj < 0.05 | Default works |
| Want more power at same FDR | `results(dds, filterFun = ihw)` | IHW typically gains 5-20% |
| Pre-specified biological fold-change matters | `results(dds, lfcThreshold = log2(1.5), altHypothesis = 'greaterAbs')` OR `glmTreat(fit, lfc = log2(1.5))` | Post-hoc `padj<0.05 & abs(LFC)>1` does NOT control FDR for the magnitude claim |
| Ranking for GSEA preranked | `stat` (Wald Z) for DESeq2 OR shrunken LFC | Never use unshrunken LFC -- low-count noise dominates |
| ORA input | Subset by padj<0.05; background = ALL TESTED genes (post-independent-filtering) | Background = "all genes in genome" is wrong; pre-filtering already excluded many |
| Many NA padj including biologically interesting genes | Diagnose: independent filtering vs Cook's vs all-zero; turn off the offending filter only for that gene set | Blanket `na.omit` discards signal |
| Multi-condition design | LRT for "any change" first; pairwise per-level Wald for effect sizes | LRT padj is omnibus; LRT LFC is one specific coefficient |
| Small n (<=3/group) | Report as exploratory, top hits only | Schurch 2016: tools miss 20-40% of true positives at n=3 |
| Human / mouse with mixed sexes | Include sex as covariate; run sex-stratified sensitivity | SABV mandate; sex effect is real and chromosomal |
| Prokaryotic | Use Prokka/Bakta GFF, KEGG strain code | Ensembl/org.db are eukaryote-only |

## Extracting Results

**Goal:** Pull DE estimates and p-values from a fitted DESeq2 or edgeR object into a usable data frame with explicit contrast naming.

**Approach:** `results()` (DESeq2) or `topTags()` (edgeR) with explicit `name=` / `coef=`; convert to data.frame; preserve row order if planning to join with annotation.

```r
library(DESeq2)
library(dplyr)

resultsNames(dds)
res <- results(dds, name = 'condition_treated_vs_control', alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
```

```r
library(edgeR)
tt <- topTags(qlf, n = Inf, sort.by = 'none')$table
tt$gene <- rownames(tt)
```

`sort.by = 'none'` in `topTags` preserves the original gene order -- critical when joining with an annotation table by row index. Default is sort by p-value.

Column name reminder (a recurring cross-tool bug):

| Tool | LFC column | Adjusted p-value column |
|------|------------|------------------------|
| DESeq2 | `log2FoldChange` | `padj` |
| edgeR | `logFC` | `FDR` |
| limma `topTable` | `logFC` | `adj.P.Val` |
| limma `topTreat` | `logFC` | `adj.P.Val` (post-TREAT) |

## TREAT vs Post-hoc LFC Filtering

**Goal:** Make a defensible FDR claim about "biologically meaningful fold change" genes.

**Approach:** Use TREAT or `lfcThreshold=` to test a magnitude hypothesis with proper FDR control. Post-hoc filtering of `padj<0.05 & abs(LFC)>tau` does NOT control FDR for the magnitude claim.

```r
res_treat <- results(dds, lfcThreshold = log2(1.5), altHypothesis = 'greaterAbs', alpha = 0.05)

# edgeR equivalent
tr <- glmTreat(fit, coef = 2, lfc = log2(1.5))
```

What a reviewer is really probing with a "200 genes >2x changed at FDR 5%" claim: is the FDR for the change>2x claim or for the change-non-zero claim? Post-hoc filtering controls FDR only for the latter. TREAT (or `lfcThreshold=`) controls FDR for the former. McCarthy & Smyth 2009 *Bioinformatics* 25:765 is the canonical citation.

## IHW for Better Power

**Goal:** Gain 5-20% more discoveries at the same FDR by weighting p-values with a covariate (typically baseMean) that informs power but is independent of the null.

**Approach:** `results(dds, filterFun = ihw)` replaces independent filtering with Ignatiadis 2016 hypothesis weighting.

```r
library(IHW)
res_ihw <- results(dds, filterFun = ihw, alpha = 0.05)
```

When IHW does NOT help:
- Small number of tests (<5000 after filtering) -- not enough data to learn weights
- Covariate is treatment-correlated (violates the null-independence requirement)
- Covariate uninformative about test power

Storey q-value as an alternative (different framework -- estimates pi_0 fraction of true nulls):

```r
library(qvalue)
qv <- qvalue(res$pvalue[!is.na(res$pvalue)])
res$qvalue <- NA
res$qvalue[!is.na(res$pvalue)] <- qv$qvalues
```

ashr lfsr (local false sign rate -- probability the estimated direction is wrong):

```r
res_ashr <- lfcShrink(dds, coef = 'condition_treated_vs_control',
                       type = 'ashr', svalue = TRUE)
res_ashr$svalue  # FDR-like, based on lfsr; requires svalue=TRUE
```

`svalue=TRUE` is required to populate the `svalue` column; the default returns the standard `pvalue`/`padj` columns only. lfsr and padj are NOT interchangeable. lfsr asks "P(sign wrong)"; padj asks "expected fraction of false discoveries". When reporting, state which.

## P-value Histogram Diagnostics

**Goal:** Diagnose model misspecification, hidden batch effects, or over-correction by inspecting the raw p-value distribution.

**Approach:** Plot raw p-values; under a correctly specified null, the histogram is uniform with an upward spike near zero (the true DE genes).

```r
library(ggplot2)
ggplot(res_df, aes(x = pvalue)) +
    geom_histogram(bins = 50, fill = 'steelblue', color = 'white') +
    labs(x = 'P-value', y = 'Frequency', title = 'P-value distribution') +
    theme_bw()
```

| Shape | Meaning | Action |
|-------|---------|--------|
| Uniform + spike near 0 | Correct: null genes uniform, true DE near 0 | Proceed |
| Anti-conservative (U-shape; both ends spiked) | Hidden batch effect, unmodeled confounder, dispersion misspecified | Inspect PCA for batch; add covariate; check `plotDispEsts` |
| Conservative (depleted near 0, spike near 1) | Over-correction; too many covariates; wrong dispersion | Simplify model; check dispersion plot for excess shrinkage |
| Spike only at p = 1 | Discrete artifact from very-low-count genes | Pre-filter more aggressively |
| Bimodal with spike at 0.5 | Unusual; suggests a discrete categorical test masquerading | Investigate |

The histogram is one of the cheapest sanity checks in a DE pipeline; always plot it before believing the gene list.

## Filtering and Ordering

**Goal:** Subset to significant genes and rank by p-value, fold change, or expression level for downstream use.

**Approach:** dplyr-style filter + arrange; handle NA padj explicitly per the three-meanings table at the top.

```r
sig <- res_df %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1, baseMean > 10) %>%
    arrange(padj)

# Up- vs down-regulated
up   <- sig %>% filter(log2FoldChange > 0)
down <- sig %>% filter(log2FoldChange < 0)

# Summary
n_tested <- sum(!is.na(res$padj))
n_sig    <- sum(res$padj < 0.05, na.rm = TRUE)
cat(sprintf('Tested: %d   Significant (padj<0.05): %d   Up: %d   Down: %d\n',
            n_tested, n_sig, sum(sig$log2FoldChange > 0), sum(sig$log2FoldChange < 0)))
```

## Gene Annotation

**Goal:** Map gene IDs to symbols, descriptions, and cross-database identifiers for human-readable results.

**Approach:** Prefer `AnnotationDbi::mapIds` with org.db (fast, local, version-pinned); fall back to biomaRt or mygene for symbols/aliases not in org.db; for prokaryotes, use Prokka/Bakta GFF.

```r
library(org.Hs.eg.db)
library(AnnotationDbi)

res_df$symbol <- mapIds(org.Hs.eg.db, keys = sub('\\..*', '', res_df$gene),
                         keytype = 'ENSEMBL', column = 'SYMBOL', multiVals = 'first')
res_df$entrez <- mapIds(org.Hs.eg.db, keys = sub('\\..*', '', res_df$gene),
                         keytype = 'ENSEMBL', column = 'ENTREZID', multiVals = 'first')
```

The `sub('\\..*', '', ...)` strips the Ensembl version. CAUTION: this regex destroys the `_PAR_Y` suffix in GENCODE 25-43 PAR genes -- use `sub('\\.[0-9]+(_PAR_Y)?$', '\\1', ...)` to preserve. See `expression-matrix/gene-id-mapping` for full details.

For HGNC symbols changed since 2020 (`SEPT1` -> `SEPTIN1`, `MARCH1` -> `MARCHF1`, `MARC1` -> `MTARC1`, `DEC1` -> `DELEC1`) old symbol-keyed downstream tools silently drop genes. Always join on stable Ensembl or Entrez IDs; use symbols as display labels only.

For prokaryotes:

```r
library(rtracklayer)
gff <- import('annotation.gff3')
gene_info <- as.data.frame(gff[gff$type == 'gene',
                                c('locus_tag', 'Name', 'product')])
res_annotated <- merge(res_df, gene_info, by.x = 'gene',
                        by.y = 'locus_tag', all.x = TRUE)
```

## GSEA Preranked Input

**Goal:** Produce a ranked list of all genes (no significance filter) for fgsea / clusterProfiler GSEA.

**Approach:** Rank by Wald statistic (DESeq2 `stat`) or shrunken LFC. NEVER use a filtered set as GSEA input -- GSEA's permutation null requires the full background.

```r
gsea_ranks <- res_df$stat
names(gsea_ranks) <- res_df$gene
gsea_ranks <- sort(gsea_ranks[!is.na(gsea_ranks)], decreasing = TRUE)

# edgeR equivalent
gsea_ranks_edger <- sign(tt$logFC) * -log10(tt$PValue)
names(gsea_ranks_edger) <- rownames(tt)
gsea_ranks_edger <- sort(gsea_ranks_edger[is.finite(gsea_ranks_edger)],
                         decreasing = TRUE)
```

`stat` (Wald Z) is preferred over raw LFC for GSEA because it combines effect and precision in one number. Unshrunken LFC is dominated by low-count noise.

## ORA Input

**Goal:** Run over-representation analysis (enrichGO, enrichKEGG) on a significant gene list with the correct background.

**Approach:** Subset to padj<0.05; background = ALL TESTED genes (post-independent-filtering), NOT the genome.

```r
library(clusterProfiler)

sig_entrez <- na.omit(res_df$entrez[res_df$padj < 0.05])
bg_entrez  <- na.omit(res_df$entrez[!is.na(res_df$padj)])

ora <- enrichGO(gene          = sig_entrez,
                universe      = bg_entrez,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = 'BP',
                pAdjustMethod = 'BH')
```

Common mistake: omitting `universe=` lets clusterProfiler default to "all annotated genes for this organism" -- which includes thousands of genes never tested. The resulting enrichment p-values are wrong (too small). The background MUST be the tested set.

## Cross-Tool Concordance Check

```r
deseq2_sig <- rownames(subset(deseq2_res, padj < 0.05))
edger_sig  <- rownames(subset(edger_tt,  FDR  < 0.05))

common      <- intersect(deseq2_sig, edger_sig)
deseq2_only <- setdiff(deseq2_sig, edger_sig)
edger_only  <- setdiff(edger_sig, deseq2_sig)

cat(sprintf('DESeq2 sig: %d   edgeR sig: %d   Common: %d (%.1f%%)\n',
            length(deseq2_sig), length(edger_sig), length(common),
            100 * length(common) / min(length(deseq2_sig), length(edger_sig))))
```

Concordance >70% at the top 500: robust. <60%: suspect filtering, normalization, or design difference -- not a tool difference. Run both pipelines with the same filtering and design to isolate.

## Per-Method Failure Modes

### Dropped a key gene by removing NAs

**Trigger:** Pipeline does `res_df <- na.omit(res_df)`; downstream gene of interest is missing from results.

**Mechanism:** Gene was flagged by independent filtering OR Cook's distance; padj is NA but the biology is real.

**Symptom:** A gene with clear differential expression in the count matrix is absent from the results table.

**Fix:** Diagnose which filter fired (independent filtering vs Cook's vs all-zero); rerun `results()` with the appropriate filter off (`independentFiltering = FALSE` or `cooksCutoff = FALSE`).

### Reported FDR on a magnitude-filtered gene set

**Trigger:** Methods section says "genes with padj < 0.05 and abs(LFC) > 1 (FDR < 5%)".

**Mechanism:** BH controls FDR for the |LFC| > 0 hypothesis, not the |LFC| > 1 hypothesis. The post-hoc filter adds no FDR control.

**Symptom:** Reviewer challenges the FDR claim; replication studies show many of the filtered genes are not the magnitude expected.

**Fix:** Use TREAT (`glmTreat`) or `lfcThreshold=` to test the magnitude hypothesis with proper FDR control. Re-do the methods sentence to match what was actually computed.

### ORA universe wrong

**Trigger:** ORA p-values look implausibly small for a small significant gene set.

**Mechanism:** `universe=` argument omitted; clusterProfiler defaulted to all annotated genes in the organism, including thousands never in the tested set.

**Symptom:** Many enriched pathways at strict thresholds; results don't replicate; reviewer questions the background.

**Fix:** Explicitly pass `universe = bg_entrez` where `bg_entrez` is the set of tested gene IDs (i.e., those with non-NA padj).

### "Significant" gene list in n=3 study doesn't replicate

**Trigger:** Small RNA-seq study finds 200 DE genes; validation in independent cohort recovers 60.

**Mechanism:** Schurch 2016 *RNA* 22:839: at n=3/group, all tools miss 20-40% of true positives compared to n=30. Variability of the gene list itself is high.

**Symptom:** 30-50% replication of the gene list across independent runs of the SAME data.

**Fix:** Frame the small-n DE list as hypothesis-generating, not as a stable set of facts. Validate top hits orthogonally before drawing conclusions. Use TREAT for biologically meaningful thresholds to require larger effects.

### Sex-confounded design gives spurious chrX/chrY signal

**Trigger:** Mixed-sex cohort; sex not in the design; many chrY genes call as DE.

**Mechanism:** Sex distribution differs across the experimental groups; the "treatment effect" partially captures sex.

**Symptom:** chrY genes (DDX3Y, RPS4Y1, UTY) and XIST dominate the top DE list.

**Fix:** Include sex in the design (`~ sex + condition`); rerun. For chrX/chrY-specific analyses, sex MUST be in the model or the analysis is uninterpretable. Mauvais-Jarvis et al. 2020 *Lancet* 396:565 reviews the SABV requirement.

## Common errors

| Error / symptom | Cause | Fix |
|-----------------|-------|-----|
| `$FDR` not found on DESeq2 result | DESeq2 uses `padj`; edgeR uses `FDR` | Check tool, use correct column |
| `summary(res)` shows different cutoff than `results(alpha=)` | `summary(res, alpha=)` defaults to 0.1 | Pass `alpha` explicitly to `summary()` |
| All `padj` NA | All genes filtered (rare; usually a data problem) | Check `independentFilteringResults(res)`; inspect baseMean distribution |
| Direction of LFC reversed | Reference level not set; alphabetical default | `relevel()` BEFORE `DESeq()` |
| Gene symbol mapping rate <50% | Mixed Ensembl versions; recent HGNC renames | Verify Ensembl release, check for SEPT/MARCH/MARC renames |
| `enrichGO` reports thousands of pathways | Wrong `universe=` | Pass `universe = bg_entrez` (tested set, not genome) |

## References

- Benjamini Y, Hochberg Y. 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc Ser B* 57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x
- Storey JD. 2003. The positive false discovery rate: a Bayesian interpretation and the q-value. *Ann Stat* 31(6):2013-2035. doi:10.1214/aos/1074290335
- Ignatiadis N, Klaus B, Zaugg JB, Huber W. 2016. Data-driven hypothesis weighting increases detection power in genome-scale multiple testing. *Nat Methods* 13(7):577-580. doi:10.1038/nmeth.3885
- Bourgon R, Gentleman R, Huber W. 2010. Independent filtering increases detection power for high-throughput experiments. *PNAS* 107(21):9546-9551. doi:10.1073/pnas.0914005107
- Stephens M. 2017. False discovery rates: a new deal. *Biostatistics* 18(2):275-294. doi:10.1093/biostatistics/kxw041
- McCarthy DJ, Smyth GK. 2009. Testing significance relative to a fold-change threshold is a TREAT. *Bioinformatics* 25(6):765-771. doi:10.1093/bioinformatics/btp053
- Schurch NJ et al. 2016. How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? *RNA* 22(6):839-851. doi:10.1261/rna.053959.115
- Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15(12):550. doi:10.1186/s13059-014-0550-8
- Mauvais-Jarvis F et al. 2020. Sex and gender: modifiers of health, disease, and medicine. *Lancet* 396(10250):565-582. doi:10.1016/S0140-6736(20)31561-0
- Bruford EA et al. 2020. Guidelines for human gene nomenclature. *Nat Genet* 52:754-758. doi:10.1038/s41588-020-0669-3
- Ziemann M, Eren Y, El-Osta A. 2016. Gene name errors are widespread in the scientific literature. *Genome Biol* 17:177. doi:10.1186/s13059-016-1044-7
- Wu T et al. 2021. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation* 2(3):100141. doi:10.1016/j.xinn.2021.100141

## Related Skills

- deseq2-basics - Generate DESeq2 results; design, contrasts, LRT, shrinkage
- edger-basics - Generate edgeR results; QL F-test, TREAT, voom
- de-visualization - P-value histogram, MA plot, volcano with shrunken LFC, heatmap, sample distance
- batch-correction - Include batch in design (vs Nygaard 2016 cardinal sin)
- timeseries-de - LRT-with-reduced-model patterns for time
- expression-matrix/gene-id-mapping - ID conversion, HGNC renames, ortholog mapping
- expression-matrix/metadata-joins - Sex covariate, paired design, sample swap detection
- pathway-analysis/go-enrichment - ORA with proper background
- pathway-analysis/gsea - GSEA preranked input from `stat` or shrunken LFC
- pathway-analysis/kegg-pathways - KEGG with strain-specific organism codes
- data-visualization/volcano-and-ma-plots - Custom volcano with apeglm-shrunken LFC
