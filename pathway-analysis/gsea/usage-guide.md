# GSEA - Usage Guide

## Overview
Gene Set Enrichment Analysis (GSEA) tests whether a gene set shifts coordinately toward one end of a genome-wide ranked gene list, without ever choosing a significance cutoff. The result is a deterministic function of three usually-implicit choices: the ranking statistic, the weight exponent p, and which labels are permuted. This skill owns the threshold-free Functional Class Scoring machinery (the running-sum enrichment score, the ranking-metric choice, the permutation null, NES/FDR, the leading edge) plus the per-sample ssGSEA/GSVA variants, using clusterProfiler's gseGO/gseKEGG/gsePathway/GSEA on the fgseaMultilevel engine.

## Prerequisites
```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db', 'fgsea'))
install.packages('msigdbr')                    # MSigDB gene sets
BiocManager::install(c('ReactomePA', 'GSVA'))  # Reactome GSEA; per-sample scores
```

Conceptual prerequisites:
- The input is a NAMED numeric vector sorted strictly DECREASING - one statistic per gene, deduplicated, NAs removed. An unranked gene list is the wrong input (route to go-enrichment).
- Rank by a signed, variance-calibrated metric: DESeq2 `stat`, limma moderated `t`, or `sign(log2FC) * -log10(p)` for edgeR. Never rank by a raw p-value (it erases direction) or bare log2FC from raw counts (low-count outliers hijack the leading edge).
- Gene IDs must match the method: ENTREZ-style for gseKEGG/gsePathway/MSigDB. The DE list and ranking statistic come from differential-expression/de-results.
- gseKEGG queries the live KEGG REST API and needs internet at run time; its results are not reproducible across KEGG releases (pin the run date). gseGO, gsePathway, and MSigDB GSEA use local annotation and are reproducible given the package version.
- Always `set.seed()` before any GSEA run; the multilevel Monte Carlo is stochastic.

## Quick Start
Tell your AI agent what you want to do:
- "Run GSEA on my DESeq2 results ranked by the Wald statistic"
- "Find GO biological processes with coordinated expression changes across all my genes"
- "Run GSEA against the MSigDB Hallmark collection on my ranked list"
- "Score each sample for pathway activity with GSVA for downstream clustering"

## Example Prompts

### Preranked GSEA on GO
> "I have a full DESeq2 result with the Wald statistic for every gene. Build a named decreasing ranked vector from the `stat` column, run GO biological-process GSEA with a fixed seed, and give me the top terms by adjusted p-value with their NES and leading-edge genes."

### ORA vs GSEA choice
> "I have a complete ranked DE result for all genes and no obvious significance cutoff. Decide whether ORA or GSEA is appropriate and run whichever fits on GO terms, explaining the choice."

### Ranking metric
> "Run GSEA using a signed p-value (`sign(log2FC) * -log10(p)`) as the ranking statistic instead of fold change, because my data came from edgeR which has no Wald-equivalent column."

### MSigDB and KEGG
> "Run GSEA against the MSigDB Hallmark collection on my ranked human gene list, then also run gseKEGG and note that KEGG queries the live database."

### Per-sample scores
> "Convert my expression matrix into a per-sample pathway-activity matrix with GSVA so I can cluster samples and correlate the scores with survival."

## What the Agent Will Do
1. Load the DE results and select a signed, variance-calibrated ranking metric matched to the DE tool.
2. Build a named numeric vector, drop NAs, deduplicate to one statistic per gene, and sort strictly decreasing.
3. Convert gene IDs to the type the chosen database expects (ENTREZ for KEGG/Reactome/MSigDB).
4. Fix the seed, set `eps=0` for exact tiny p-values, and run gseGO/gseKEGG/gsePathway or GSEA with a TERM2GENE collection.
5. Report `p.adjust`/`qvalue` and NES, then read `core_enrichment` (the leading edge) as the interpretable core.
6. For per-sample activity, build a GSVA/ssGSEA parameter object and produce a gene-set-by-sample score matrix.

## GSEA vs Over-Representation

| Feature | Over-representation (ORA) | GSEA |
|---------|--------------------------|------|
| Input | pre-selected gene list + a background | ranked vector of ALL genes |
| Cutoff | requires a significance threshold | no arbitrary cutoff |
| Detects | strong individual changes | coordinated subtle shifts |
| Functions | enrichGO, enrichKEGG | gseGO, gseKEGG, GSEA |

The full ORA-vs-FCS-vs-topology fork and the competitive/self-contained null theory live in the category README.

## Choosing a Ranking Statistic

| DE source | Ranking metric | Column | Notes |
|-----------|----------------|--------|-------|
| DESeq2 | Wald statistic | `stat` | best single choice for RNA-seq; signed + variance-calibrated |
| limma / voom | moderated t-statistic | `t` | empirical-Bayes shrinkage; signed |
| edgeR | signed p-value | `sign(logFC) * -log10(PValue)` | no Wald-equivalent; clamp p==0 with pmax(p, 1e-300) |
| any tool, last resort | log2FC alone | `log2FoldChange` | magnitude only; noisy for low-count genes |

apeglm/ashr-shrunk DESeq2 results drop the `stat` column - pull `stat` from the unshrunk `results(dds)` if ranking by it. Never use `lfcShrink(type='normal')` for ranking.

## Interpreting NES (Normalized Enrichment Score)
- Positive NES: the set piles up at the top of the ranking (up in the contrast).
- Negative NES: the set piles up at the bottom (down in the contrast).
- Check FDR (`p.adjust`) first; a high |NES| with non-significant FDR is meaningless.
- clusterProfiler `p.adjust` is a BH-adjusted FDR, computed differently from the Broad desktop tool's empirical-null FDR - the "FDR < 0.25" convention belongs to the latter, so a defensible BH cutoff is usually 0.05; reserve 0.25 for genuinely exploratory work and say which was used.
- Read the leading edge (`core_enrichment`): a large, concentrated leading edge means coordinated regulation; a 1-2 gene leading edge means outliers, not a pathway.

## Tips
- The ranking is the experiment; the gene sets are just the question asked of it. A bad ranking is faithfully reported as a ranking artifact.
- Preranked GSEA (clusterProfiler/fgsea) uses GENE permutation, which destroys gene-gene correlation and is anti-conservative for co-regulated sets. Report the permutation type; when the design matrix is available, use CAMERA (`limma::camera`) for a correlation-honest competitive test.
- `nPerm` no longer exists; tiny-p accuracy is governed by `eps` (set `eps=0` for exact). The engine is `fgseaMultilevel` via `by='fgsea'` (the default).
- ssGSEA and GSVA are NOT a contrast test - they produce a per-sample activity matrix with no per-set p-value. GSVA >= 1.50 needs the parameter-object API: `gsva(gsvaParam(expr, sets))`.
- msigdbr 26.x renamed `category=` to `collection=` and `gs_cat` to `gs_collection`; the Entrez column is `ncbi_gene` (older releases used `entrez_gene`). Check `?msigdbr` and `names()` for the installed version.
- If no terms are enriched, check the ranking metric (is it signed?), confirm the vector is named and sorted decreasing, deduplicate gene IDs, and verify the ID type matches the database.
- See enrichment-visualization for gseaplot2(), ridgeplot(), and dotplot() of GSEA results.

## Related Skills

- go-enrichment - Gene-list ORA alternative when no ranking exists
- kegg-pathways - KEGG pathway and module enrichment and ID conventions
- reactome-pathways - Reactome curated-pathway enrichment (local DB)
- wikipathways - WikiPathways community-pathway enrichment
- enrichment-visualization - gseaplot2, ridgeplot, and dotplot of GSEA results
- differential-expression/de-results - Source of the ranking statistic and the padj column conventions
- workflows/expression-to-pathways - End-to-end DE-to-enrichment pipeline
