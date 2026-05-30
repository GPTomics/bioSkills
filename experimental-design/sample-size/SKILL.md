---
name: bio-experimental-design-sample-size
description: Estimates the minimum biological replicates (or cells/events) for a target power at a target FDR in genomics experiments using ssizeRNA, PROPER, powsimR for scRNA-seq, and pilot-data dispersion estimation from DESeq2/edgeR. Covers the biological-versus-technical replication distinction (technical replicates do not add degrees of freedom for biological inference), replicate-number-versus-sequencing-depth budgeting, scRNA-seq sample-versus-cell allocation under a pseudobulk model, and the critique that "n=3" is a publication convention rather than a power calculation. Use when budgeting a sequencing experiment, writing the sample-size justification in a grant, estimating replicates from pilot data, allocating a fixed budget between samples and depth, or planning scRNA-seq cohort size. For clinical-trial sample size see clinical-biostatistics/power-and-sample-size; for the power-given-n direction see experimental-design/power-analysis.
tool_type: r
primary_tool: ssizeRNA
---

## Version Compatibility

Reference examples tested with: ssizeRNA 1.3+, PROPER 1.34+, powsimR 1.2+ (GitHub), DESeq2 1.42+, edgeR 4.0+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws an error, introspect the installed package and adapt to the actual API. Notes: ssizeRNA provides `ssizeRNA_single()` (one mean/dispersion for all genes), `ssizeRNA_vary()` (genes vary), and `check.power()` (average power and true FDR for a given n); powsimR is GitHub-only with drifting signatures. Verify against the installed help before use.

# Sample Size for Genomics Experiments

**"How many samples do I need?"** -> Find the smallest number of biological replicates per group that achieves a target marginal power at a target FDR, given the dispersion and effect-size distribution expected for the assay — counting biological units, not measurements.
- R: `ssizeRNA::ssizeRNA_vary()`, `ssizeRNA::check.power()` — FDR-aware NB sample size; pilot dispersions from `DESeq2`/`edgeR`

## The Single Most Important Modern Insight -- The Biological Replicate Is the Unit, and n=3 Is a Convention

Sample size is a count of **biological replicates** — independent experimental units (animals, donors, cultures from independent passages), not measurements. Technical replicates (one library split across lanes, one RNA split into preps) reduce measurement noise but add **no degrees of freedom** for biological inference; averaging them into their biological unit is correct, and selling "n = 3 samples x 3 technical reps = 9" as biological power is a standard error (Blainey, Krzywinski & Altman 2014 *Nat Methods* 11:879). The ubiquitous **"n=3" is a publication convention, not a calculation**: in the 48-vs-48 yeast benchmark, **>=6** biological replicates were needed to recover most true DE genes at realistic effect sizes, and below that the choice of DE tool mattered more than at higher n (Schurch 2016 *RNA* 22:839). Human and primary material, with higher dispersion, need more. For single-cell, the corollary is sharp: population-level DE power is set by the **number of donors**, not the number of cells, because cells are pseudoreplicates — pseudobulk per donor is the correct unit (Squair 2021 *Nat Commun* 12:5692; Murphy & Skene 2022 *Nat Commun* 13:7851).

## Algorithmic Taxonomy

| Approach | Model | Tool | Strength | Fails / costs when |
|----------|-------|------|----------|--------------------|
| FDR-aware NB sample size | NB, varying mean/dispersion | `ssizeRNA::ssizeRNA_vary` | controls average power at a true FDR | needs a dispersion/expression model |
| Pilot-dispersion simulation | empirical dispersions from pilot | `PROPER`, `powsimR` | most defensible; study-specific | requires a pilot dataset |
| Single-parameter NB | one mean/dispersion for all genes | `ssizeRNA::ssizeRNA_single` | quick; transparent | ignores the mean-dispersion trend |
| Verify a planned n | average power + true FDR at fixed n | `ssizeRNA::check.power` | sanity-checks a budget-driven n | not a search over n |
| scRNA-seq cohort sizing | pseudobulk over donors | `powsimR` | counts the right unit (donors) | cell-level sizing is wrong unit |
| Per-feature t-test n | Gaussian (Cohen's d) | `pwr::pwr.t.test` | proteomics/continuous after transform | wrong for raw counts |

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Bulk RNA-seq, pilot available | estimate dispersions, then `ssizeRNA_vary`/PROPER | study-specific dispersion beats a guess |
| Bulk RNA-seq, no pilot | `ssizeRNA_vary` with a literature dispersion, stated as approximate | transparent starting point |
| Budget already fixed at some n | `check.power` to report achieved power and true FDR | answers "is this n adequate?" |
| scRNA-seq disease vs control | size the number of DONORS (pseudobulk; powsimR) | population power scales with donors |
| ChIP/ATAC/methylation | NB sample size per region; assay floor as minimum | overdispersed counts; detection floor |
| Proteomics (continuous) | `pwr::pwr.t.test` per protein, with missingness caveat | Gaussian after transform |
| Have technical replicates | collapse to biological units first | technical reps add no biological df |
| Clinical-trial endpoint | -> clinical-biostatistics/power-and-sample-size | regulated regime |

## FDR-Aware NB Sample Size -- ssizeRNA

**Goal:** Find the minimum biological replicates per group for a target power at a target FDR, accounting for the proportion of DE genes and the mean-dispersion structure.

**Approach:** Specify the number of genes, the proportion non-DE (pi0), the mean count and dispersion (ideally from pilot data), the fold change, the target FDR, and the target power; let `ssizeRNA_vary` search replicate numbers and return the smallest that reaches the target.

```r
library(ssizeRNA)
res <- ssizeRNA_vary(nGenes = 20000, pi0 = 0.95,        # 5% DE
                     mu = 10, disp = 0.2,                # mean count + dispersion (from pilot ideally)
                     fc = 1.5, fdr = 0.05, power = 0.80,
                     maxN = 30)
res$ssize                                                # minimum n per group

# Verify a budget-fixed n: average power and TRUE realized FDR
check.power(nGenes = 20000, pi0 = 0.95, m = 6, mu = 10, disp = 0.2, fc = 1.5, fdr = 0.05, sims = 50)
```

## Pilot Dispersions Drive Honest Sample Size

**Goal:** Replace a guessed CV with a measured dispersion-mean trend from pilot data.

**Approach:** Fit dispersions on the pilot with DESeq2 or edgeR, summarize them, and feed them into the simulation-based estimator (PROPER or powsimR) rather than a single-CV closed form.

```r
library(DESeq2)
dds <- DESeqDataSetFromMatrix(pilot_counts, pilot_coldata, ~ condition)
dds <- DESeq(dds)
disp <- dispersions(dds)                                 # per-gene dispersion estimates
summary(disp[is.finite(disp)])                           # feed median/trend to PROPER/powsimR
# A literature CV can be off by ~2x; a pilot dispersion is the defensible input.
```

## Biological vs Technical Replication

Technical replicates estimate measurement variance; biological replicates estimate the variance that generalizes to the population, and only the latter supports inference about the biology. Average or sum technical replicates into their biological unit before any test. "n = 3 samples x 3 technical reps" is n = 3, not n = 9 (Blainey 2014). This is the sample-size face of the experimental-unit principle (see experimental-design/randomization-blocking).

## Replicates vs Depth Under a Fixed Budget

Once depth is adequate (roughly >=10-20M mapped reads for bulk RNA-seq DE), additional biological replicates buy more power than additional depth (Liu 2014 *Bioinformatics* 30:301). Allocate a fixed budget toward more biological units first. scRNA-seq has an analogous rule at the donor level: more donors beat more cells per donor for population DE, with cells per cell type showing diminishing returns past a few hundred (Squair 2021; Murphy-Skene 2022).

## Sample Size by Assay (floors under favorable conditions, not targets)

| Assay | Practical minimum | For small effects | Source / note |
|-------|-------------------|-------------------|---------------|
| Bulk RNA-seq | 3 (convention) | 6-12 | Schurch 2016 *RNA* 22:839: >=6 recovers most true DE |
| scRNA-seq (population DE) | 3 donors | 6+ donors | Squair 2021; donors, not cells, drive power |
| ATAC-seq | 2 | 4-6 | library complexity + peak detection floor |
| ChIP-seq | 2 | 3-4 | IDR reproducibility framework (ENCODE) |
| Proteomics (DIA/TMT) | 3 | 6-10 | higher missingness; MNAR |
| Methylation (array/WGBS) | 4 | 8-12 | high per-CpG variance |

The "minimum" columns are floors that assume low dispersion and large effects; treat them as the smallest defensible n only after a pilot or literature dispersion supports them.

## Per-Method Failure Modes

### Technical replicates counted as biological n
- **Trigger:** "n = 9: 3 samples x 3 technical reps."
- **Mechanism:** technical reps add no biological degrees of freedom (Blainey 2014).
- **Symptom:** over-stated power; results do not generalize.
- **Fix:** collapse technical reps to the biological unit; biological n = 3.

### n=3 by convention
- **Trigger:** choosing 3 because "everyone uses 3."
- **Mechanism:** 3 is a habit, not a calculation; misses many true DE (Schurch 2016).
- **Symptom:** chronic underpowering, irreproducibility.
- **Fix:** size from dispersion + target FDR; expect >=6 for realistic effects, more for human material.

### scRNA-seq sized on cells
- **Trigger:** "100k cells from 2 donors is plenty."
- **Mechanism:** population power scales with donors; cells are pseudoreplicates (Squair 2021).
- **Symptom:** false-discovery-laden DE that does not replicate.
- **Fix:** budget for more donors; size on a pseudobulk model.

### Guessed CV instead of pilot dispersion
- **Trigger:** "human samples are ~0.4, so use 0.4."
- **Mechanism:** real dispersion is study-specific; the guess can be off ~2x.
- **Symptom:** the planned n is wrong by a large factor.
- **Fix:** estimate dispersion from any available pilot (DESeq2/edgeR).

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| >=6 biological replicates for bulk RNA-seq DE | Schurch 2016 *RNA* 22:839 | recovers most true DE at realistic effects |
| n=3 is a convention, not a calculation | Schurch 2016 | low power and tool-dependent below 6 |
| Donors, not cells, set scRNA-seq DE power | Squair 2021 *Nat Commun* 12:5692 | cells are pseudoreplicates |
| Technical reps add 0 biological df | Blainey 2014 *Nat Methods* 11:879 | only biological reps generalize |
| Depth saturates ~10-20M reads; add replicates | Liu 2014 *Bioinformatics* 30:301 | biological variance dominates |
| Add 10-20% extra units for failures | common practice | RNA degradation, failed libraries |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Over-stated power | technical reps counted as n | collapse to biological units |
| Underpowered at n=3 | convention not calculation | size to >=6 (or pilot-driven) |
| scRNA-seq DE does not replicate | sized on cells | size on donors (pseudobulk) |
| Planned n off by a large factor | guessed CV | estimate dispersion from pilot |
| Study fails after sample loss | no failure margin | add 10-20% extra units |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Why this n?" | smallest n reaching marginal power >= 0.8 at FDR 0.05 for the minimum meaningful FC; power curve provided |
| "Where did dispersion come from?" | estimated from pilot (DESeq2); literature value used only as a cross-check |
| "Is n=3 enough?" | no; sized to >=6 per Schurch 2016 for realistic effects |
| "Why so many donors for scRNA-seq?" | population DE power scales with donors, not cells (Squair 2021) |
| "Technical replicates?" | collapsed to biological units; they add no biological degrees of freedom |

## References

- Bi R, Liu P. 2016. Sample size calculation while controlling false discovery rate for differential expression analysis with RNA-sequencing experiments. *BMC Bioinformatics* 17:146.
- Schurch NJ, Schofield P, Gierliński M, et al. 2016. How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? *RNA* 22:839-851.
- Blainey P, Krzywinski M, Altman N. 2014. Points of significance: replication. *Nat Methods* 11:879-880.
- Liu Y, Zhou J, White KP. 2014. RNA-seq differential expression studies: more sequence or more replication? *Bioinformatics* 30:301-304.
- Squair JW, Gautier M, Kathe C, et al. 2021. Confronting false discoveries in single-cell differential expression. *Nat Commun* 12:5692.
- Murphy AE, Skene NG. 2022. A balanced measure shows superior performance of pseudobulk methods in single-cell RNA-sequencing analysis. *Nat Commun* 13:7851.
- Wu H, Wang C, Wu Z. 2015. PROPER: comprehensive power evaluation for differential expression using RNA-seq. *Bioinformatics* 31:233-241.
- Vieth B, Ziegenhain C, Parekh S, Enard W, Hellmann I. 2017. powsimR: power analysis for bulk and single cell RNA-seq experiments. *Bioinformatics* 33:3486-3488.

## Related Skills

- power-analysis - The power-given-n direction and simulation-based power
- randomization-blocking - The experimental unit defines what is counted as a replicate
- batch-design - Balanced designs assume equal n per group
- differential-expression/deseq2-basics - Estimating pilot dispersions for the sample-size model
- single-cell/preprocessing - Pseudobulk aggregation underlying scRNA-seq cohort sizing
- clinical-biostatistics/power-and-sample-size - Sample size for regulated clinical trials
