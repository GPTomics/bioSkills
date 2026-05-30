---
name: bio-flow-cytometry-differential-analysis
description: Differential abundance (DA) and differential state (DS) analysis for flow and mass cytometry - tests which cell populations change in frequency or marker expression between conditions using diffcyt (edgeR/voom/GLMM for DA, limma/LMM for DS), with cydar, CITRUS, and compositional methods (sccomp, scCODA, DCATS) as alternatives. Covers the sample-is-the-experimental-unit principle, design/contrast and mixed-model formulas, compositionality of cluster proportions, and FDR across clusters. Use when comparing populations between groups, choosing a DA method, handling paired/batch designs, or deciding whether compositional correction is needed.
tool_type: r
primary_tool: diffcyt
---

## Version Compatibility

Reference examples tested with: diffcyt 1.22+, CATALYST 1.26+, edgeR 4.0+, limma 3.58+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

`testDA_edgeR`/`testDS_limma` are diffcyt functions operating on count/median objects from `calcCounts`/`calcMedians`; the CATALYST-integrated path is the `diffcyt()` wrapper on the SCE. Confirm the signature with `?diffcyt` before relying on it.

# Differential Analysis

**"Compare cell populations between my conditions"** -> Test cluster frequencies (DA) and within-cluster marker expression (DS) between groups, with the sample (not the cell) as the unit.
- R: `diffcyt::diffcyt(sce, analysis_type='DA', method_DA='diffcyt-DA-edgeR', design, contrast)`
- R: `diffcyt(sce, analysis_type='DS', method_DS='diffcyt-DS-limma', ...)`

## The Single Most Important Modern Insight -- The Sample Is the Experimental Unit, Not the Cell

Tens of thousands of cells from one donor are technical PSEUDOREPLICATES, not independent observations. A per-cell test (Wilcoxon across all cells) treats them as n = cells and produces astronomically significant p-values from two mice - it is the single most common statistical sin in modern cytometry (Hurlbert 1984 *Ecol Monogr* 54:187; the cytometry mirror of the scRNA-seq pseudobulk lesson). The correct unit is the SAMPLE/subject: diffcyt aggregates cells to PER-SAMPLE-PER-CLUSTER counts (DA) and PER-SAMPLE-PER-CLUSTER arcsinh-MEDIANS (DS), then tests across samples with edgeR/limma/GLMM (Weber 2019 *Commun Biol* 2:183). Biological replication is mandatory (>= 2-3 per group); DA from a single sample per condition has no valid test. Paired with this: cluster proportions are COMPOSITIONAL (they sum to 1), so a real increase in one population mechanically forces apparent depletion in others - a source of false DA in "unchanged" clusters.

## DA vs DS, and the type/state marker link

- **DA** (differential abundance): does a cluster's FREQUENCY differ? Clusters are defined by TYPE markers.
- **DS** (differential state): within a fixed-identity cluster, does a STATE marker's expression differ? State markers were withheld from clustering for exactly this test.

## Method Taxonomy

| Method | Citation | Mechanism | When to use |
|--------|----------|-----------|-------------|
| diffcyt-DA-edgeR / voom | Weber 2019 *Commun Biol* 2:183 | edgeR/voom empirical-Bayes on per-sample counts; optional TMM | standard 2+ group with replicates (DEFAULT) |
| diffcyt-DA-GLMM / DS-LMM | Weber 2019 | random effects in the formula | paired/repeated-measures/nested (subject random effect) |
| cydar | Lun 2017 *Nat Methods* 14:707 | overlapping hyperspheres + edgeR + spatial FDR | continuum, avoid hard clusters |
| CITRUS | Bruggner 2014 *PNAS* 111:E2770 | hierarchical clustering + LASSO | predictive signature, LARGE n; correlated-not-causal; largely superseded |
| sccomp / scCODA / DCATS | Mangiola 2023 *PNAS* 120:e2203828120 / Buttner 2021 *Nat Commun* 12:6876 / Lin 2023 *Genome Biol* 24:151 | simplex-aware compositional models | strong compositional shift (one pop dominates); DCATS for assignment uncertainty |

## Run diffcyt DA and DS

**Goal:** Test abundance and state on a CATALYST-clustered SCE.

**Approach:** Build design + contrast from `ei(sce)`; the `diffcyt()` wrapper uses the stored clustering. State markers are tested in DS, type markers define DA clusters.

```r
library(CATALYST); library(diffcyt)

sce <- readRDS('sce_clustered.rds')
design   <- createDesignMatrix(ei(sce), cols_design = 'condition')
contrast <- createContrast(c(0, 1))                    # Treatment vs Control

res_DA <- diffcyt(sce, clustering_to_use = 'meta20',
                  analysis_type = 'DA', method_DA = 'diffcyt-DA-edgeR',
                  design = design, contrast = contrast)
res_DS <- diffcyt(sce, clustering_to_use = 'meta20',
                  analysis_type = 'DS', method_DS = 'diffcyt-DS-limma',
                  design = design, contrast = contrast)

library(SummarizedExperiment)
rowData(res_DA$res)        # cluster_id, logFC, p_val, p_adj (BH across clusters)
```

## Paired / Repeated-Measures (mixed models)

**Goal:** Account for within-subject correlation (e.g. pre/post on the same donor).

**Approach:** Use a GLMM/LMM method with a random effect for subject via a formula.

```r
formula <- createFormula(ei(sce), cols_fixed = 'condition', cols_random = 'patient_id')
res_DA  <- diffcyt(sce, clustering_to_use = 'meta20',
                   analysis_type = 'DA', method_DA = 'diffcyt-DA-GLMM',
                   formula = formula, contrast = createContrast(c(0, 1)))
```

## Compositional Re-Check

**Goal:** Confirm a headline single-population shift is not inducing artifactual reciprocal depletion.

**Approach:** Re-test with a simplex-aware model when one cluster changes a lot or total yield differs by group.

```r
# If a dominant population expands, the apparent depletion of others may be a simplex artifact.
# Re-test with sccomp / scCODA (reference cell type) / DCATS (assignment uncertainty)
# before reporting reciprocal depletion as independent biology.
```

## Per-Method Failure Modes

### Per-cell pseudoreplication
**Trigger:** Wilcoxon/t-test across all cells. **Mechanism:** cells aren't independent. **Symptom:** p ~ 1e-40 from few subjects. **Fix:** aggregate to per-sample summaries (diffcyt).

### Compositional false DA
**Trigger:** one population expands strongly. **Mechanism:** proportions sum to 1. **Symptom:** significant "depletion" of unrelated clusters. **Fix:** TMM only when total cell abundance is NOT itself the biological signal (else it removes real signal), or a compositional method (sccomp/scCODA/DCATS); report total-yield differences.

### Batch cleaned instead of modeled
**Trigger:** normalizing batch out then testing naively. **Mechanism:** over-correction removes real signal. **Symptom:** attenuated effects. **Fix:** include batch in the design; if batch == condition, no rescue - design it out.

### No replicates
**Trigger:** 1 sample per condition. **Mechanism:** no error term. **Symptom:** uninterpretable p. **Fix:** require >= 2-3 biological replicates per group.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| >= 2-3 biological replicates per group | Weber 2019 | minimum for a valid DA/DS error term |
| BH FDR across clusters (and clusters x markers for DS) | diffcyt | high-resolution grids have many tests |
| arcsinh median as DS statistic | Nowicka 2017 | robust per-cluster per-sample summary |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `testDA_edgeR(sce, ...)` fails | wrong signature | use the `diffcyt()` wrapper on the SCE, or `calcCounts` first |
| results empty | wrong `clustering_to_use` name | match the stored clustering id (e.g. `meta20`) |
| no DS results | state markers not flagged | set `marker_class='state'` in the panel |
| paired design ignored | used fixed-effect method | use `diffcyt-DA-GLMM` with a random effect |

## References

- Weber 2019 *Commun Biol* 2:183 — diffcyt (DA + DS).
- Bruggner 2014 *PNAS* 111(26):E2770-E2777 — CITRUS.
- Lun 2017 *Nat Methods* 14(7):707-709 — cydar hypersphere DA.
- Mangiola 2023 *PNAS* 120(33):e2203828120 — sccomp compositional analysis.
- Buttner 2021 *Nat Commun* 12:6876 — scCODA.
- Lin 2023 *Genome Biol* 24:151 — DCATS (assignment-uncertainty-aware).
- Nowicka 2017 *F1000Research* 6:748 — CyTOF workflow; arcsinh-median DS statistic.
- Hurlbert 1984 *Ecol Monogr* 54(2):187-211 — pseudoreplication.

## Related Skills

- clustering-phenotyping - Cluster (type markers) before testing
- gating-analysis - Compare manually gated population frequencies
- differential-expression/de-results - Shared edgeR/limma output semantics (padj)
- differential-expression/edger-basics - The count-model engine diffcyt reuses
- experimental-design/multiple-testing - FDR across clusters and clusters x markers
- experimental-design/batch-design - Model batch in the design, don't clean it out
