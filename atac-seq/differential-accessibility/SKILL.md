---
name: bio-atac-seq-differential-accessibility
description: Identify differentially accessible chromatin regions across conditions using DiffBind, csaw, DESeq2, or edgeR. Use when comparing ATAC-seq accessibility between treatment groups, choosing between consensus-peak vs sliding-window approaches, picking the correct normalization (full library vs reads-in-peaks), correcting batch with SVA/RUVseq, or interpreting log2FC and FDR thresholds in a chromatin context.
tool_type: r
primary_tool: DiffBind
---

## Version Compatibility

Reference examples tested with: DiffBind 3.12+, DESeq2 1.42+, edgeR 4.0+, csaw 1.36+, limma 3.58+, GenomicRanges 1.54+, ChIPseeker 1.38+, Subread 2.0+ (featureCounts), sva 3.50+, RUVSeq 1.36+.

Before using code patterns, verify installed versions match:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws unexpected errors, introspect the installed package and adapt rather than retrying.

# Differential Accessibility

**"Find chromatin regions that change accessibility between my conditions"** -> Build a sample-by-region count matrix, normalize for library size and chromatin compaction, fit a generalized linear model (negative-binomial), and extract regions with significant accessibility change.

- R (consensus-peak workflow): `DiffBind` -> count -> normalize -> contrast -> analyze
- R (window-based, no peak set): `csaw::windowCounts` + `filterWindowsGlobal` + edgeR QL F-test
- R (existing peak-count matrix): `DESeq2` or `edgeR` directly on `featureCounts` output

DiffBind is a wrapper around DESeq2 / edgeR with ATAC-aware defaults. csaw is the only peak-free option; it tests fixed-width sliding windows. The choice depends on whether peaks are stable across conditions (use DiffBind) or whether some condition has dramatically different peak structure (use csaw or rebuild consensus peaks).

## Algorithmic Taxonomy

| Tool | Model | Input | Min reps | Strength | Fails when |
|------|-------|-------|----------|----------|------------|
| DiffBind 3.x (default DESeq2) | NB GLM via DESeq2 on consensus peaks | BAM + peak files | 2-3 per group | ATAC-aware defaults; native normalization (RiP); built-in QC; blocking factors | Peaks differ dramatically between conditions (closed -> open shifts width); fewer than 2 reps per group |
| DiffBind with edgeR backend | NB GLM via edgeR-QL on consensus peaks | Same | 2-3 per group | Robust at low replicates (n=2 OK); QL test calibrates dispersion better than DESeq2 at small n | Default normalization can be off for ATAC; specify `normalize=DBA_NORM_NATIVE` |
| DESeq2 directly on peak counts | NB GLM with shrinkage | featureCounts SAF | 3+ | Maximum control; integrates with apeglm shrinkage; modern interface | Need to manually build consensus peakset; per-region pre-filter required (low counts inflate dispersion) |
| edgeR QL F-test on peak counts | NB QL (quasi-likelihood) | featureCounts | 2 | Calibrated FDR at low n (n=2 viable); robust to outlier reps | Manual consensus peakset; small library bias unless normalization explicit |
| csaw (windows) | edgeR-QL on sliding windows | BAM only | 2 | No peak set required; detects diffuse changes peaks miss; merges adjacent windows | Computationally heavy; window size choice biases results; harder to annotate downstream |
| limma-voom | linear model with mean-variance trend | log2(CPM+offset) | 3 | Fast; good calibration at moderate count | Mis-calibrated at very low counts (atac peaks often have dropouts); needs explicit voom normalization |

Methodology evolves; verify the current consensus practice (Schep & Greenleaf 2017; Reske 2020 benchmark) before locking pipelines.

## Decision Tree by Experimental Scenario

| Scenario | Recommended workflow | Why |
|----------|---------------------|-----|
| 3+ reps, similar peak structure across conditions | DiffBind (DESeq2 backend), `summits=250`, `normalize=NATIVE` | Standard pattern; peak-level inference is interpretable |
| 2 reps per condition | DiffBind with edgeR backend OR raw edgeR QL | DESeq2 underpowered at n=2; QL is robust |
| Peak structure differs dramatically (e.g., differentiation, KO of pioneer TF) | csaw windows OR rebuild consensus peakset post-hoc per condition then take union | Stable consensus peakset is invalid when chromatin landscape shifts |
| Multi-factor design (batch, sex, time) | DiffBind with `dba.contrast(..., design='~Batch + Condition')` | Standard linear model adjustment |
| Hidden batch / unknown variance | DESeq2 + SVA OR RUVseq before fitting | Empirical surrogate variables capture unknown nuisance |
| Long timecourse (5+ time points) | DESeq2 LRT (likelihood ratio test) on `~time + condition + time:condition` | Captures temporal interaction; use `differential-expression/timeseries-de` patterns |
| Diffuse / broad accessibility change (super-enhancers) | csaw with merged windows OR call broad peaks first | Narrow peaks fragment broad domains -> inflated peak count, deflated effect |
| Single-cell ATAC pseudobulk | DESeq2 on aggregated counts OR Signac::FindMarkers | See atac-seq/single-cell-atac |
| Allele-specific accessibility | csaw on heterozygous SNPs OR HOMER tagDir | Peak-level invalid because alleles share peaks |
| Plant / non-model organism | DiffBind works; just provide custom `genome` and disable annotation | Annotation step assumes UCSC TxDb; bypass if absent |

## Consensus Peak Set Strategy

The consensus peakset choice drives FDR calibration. DiffBind defaults rarely match what a chromatin biologist wants.

| Strategy | Implementation | When to use |
|----------|---------------|-------------|
| Intersection (peak in all reps) | `dba.count(minOverlap=N)` with N = total reps | Strict; for high-confidence reproducible analysis (matches IDR philosophy) |
| Union (peak in any rep) | `minOverlap=1` | Maximum sensitivity; risks single-rep artefact peaks |
| Majority rule (peak in >= half reps) | `minOverlap=ceiling(N/2)` | DiffBind default-ish; balance |
| Per-condition union, then union of unions | Compute consensus per group, then merge | Best when conditions have very different peak counts |
| Iterative overlap removal (Corces 2018) | Sort peaks by significance; greedily keep non-overlapping; fixed-width 501 bp | Standard for fixed-width consensus; required for peak-count matrices used in machine learning |

Refer to atac-seq/consensus-peakset for full coverage of fixed-width re-centering and the iterative overlap algorithm. For DiffBind, the key parameter is `summits=250` (re-center peaks on summit +/- 250 bp = 501 bp fixed width).

## Normalization: The ATAC-Specific Choice

| Normalization | DiffBind flag | What it normalizes to | When to use |
|---------------|---------------|------------------------|-------------|
| DiffBind native (RiP) | `DBA_NORM_NATIVE` | Reads in consensus peaks | Default for ATAC; protects against background drift |
| Full library size | `DBA_NORM_LIB` | Total mapped reads | When global accessibility shifts must remain visible (e.g., chromatin compaction) |
| TMM (edgeR) | `DBA_NORM_TMM` | Trimmed mean of M-values | Robust to a few highly-DA peaks dominating |
| RLE (DESeq2) | `DBA_NORM_RLE` | DESeq2 size factors | DESeq2 default; geometric-mean reference |
| Spike-in / external | none built-in; pre-scale counts | Exogenous reference (e.g., human spike in mouse) | Required when global scaling is biological |

**Trigger:** Treatment causes global chromatin compaction (e.g., HDAC inhibitor, DNMT inhibitor).

**Mechanism:** RiP normalization assumes ratio of in-peak to total reads is stable. When a treatment globally opens or closes chromatin, RiP normalizes that biology away, producing a near-null differential set.

**Symptom:** Volcano plot is symmetric about zero; PCA shows treatment effect that vanishes after normalization.

**Fix:** Use library-size normalization (`DBA_NORM_LIB`) OR use spike-in normalization (add exogenous chromatin pre-Tn5; scale by spike-in reads). Reske 2020 documented this with HDAC inhibitor.

## Per-Tool Failure Modes

### DiffBind -- Default normalization confounds global change

**Trigger:** Treatment causes whole-genome accessibility shift; cell-cycle synchronized samples; differentiation timecourse.

**Mechanism:** DiffBind default uses RiP-based normalization (`DBA_NORM_NATIVE`), which assumes total reads-in-peaks is comparable across samples. Global shift breaks this.

**Symptom:** Conditions clearly differ in PCA before normalization; after normalization PC1 is nearly noise.

**Fix:** Set `normalize=DBA_NORM_LIB` and re-run `dba.contrast` and `dba.analyze`. Re-inspect PCA; if treatment now drives PC1, the global-shift biology is preserved.

### DiffBind summits parameter -- Width-driven differential

**Trigger:** Per-rep peaks have very different widths; consensus uses union.

**Mechanism:** Without `summits=250`, DiffBind counts reads in the original peak intervals. A peak called as 200 bp in one rep and 800 bp in another inflates the count for the wider rep.

**Symptom:** Top differential peaks track peak width, not signal intensity.

**Fix:** Always set `summits=250` (or 100, depending on resolution). This re-centers all peaks on the summit and uses identical 501 (or 201) bp windows.

### csaw -- Window size and filter choice dominates results

**Trigger:** Default `width=10` (bp) windows; default `filter.global` cutoff.

**Mechanism:** Narrow windows have very low counts and inflated dispersion; default global filter discards too many windows. Results are extremely sensitive to these.

**Symptom:** Number of significant windows ranges from 200 to 200,000 across reasonable parameter sweeps.

**Fix:** Use `width=150` for ATAC (matches typical NFR fragment); set `filter.global(stat, log2(min.fold))=log2(3)` to discard low-signal windows. Validate by running on technical replicates -- ~0 differential windows is the expected outcome.

### DESeq2 -- Apeglm shrinkage with too few reps

**Trigger:** n=2 per condition; using `lfcShrink(type='apeglm')`.

**Mechanism:** Apeglm shrinks log2FC toward zero based on dispersion estimate; at n=2 dispersion is unreliable, shrinkage is over-aggressive, and biology is masked.

**Symptom:** All log2FC values cluster near zero post-shrinkage; FDR list has high p-values across the board.

**Fix:** Skip shrinkage at n=2 OR switch to edgeR QL test. If n=2 is unavoidable, report unshrunken log2FC alongside FDR; do not use shrunken FCs as the effect size.

### edgeR QL -- Filter must be aggressive enough

**Trigger:** Including peaks with mean count < 5 across all samples.

**Mechanism:** The QL F-test calibrates dispersion across all features. Including very-low-count peaks pulls dispersion estimates and inflates FDR.

**Fix:** `filterByExpr(y, group=group)` removes low-count peaks; restore peaks one at a time only if they are biologically critical and supported by at least one rep at depth.

## Reconciliation: When Tools Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| DiffBind + DESeq2 differ wildly | Different normalization (DiffBind default = RiP, DESeq2 default = RLE) | Force same normalization; differences should shrink |
| DiffBind + csaw differ | csaw catches diffuse changes peaks miss; DiffBind catches narrow peaks csaw smooths | Both can be correct; report intersection as high-confidence |
| Top hits in DiffBind have FDR > 0.5 in DESeq2 | DiffBind's blacklist filter or width re-centering changes the per-region count | Re-run DESeq2 on the exact DiffBind consensus matrix (`dba.peakset` extract) |
| Effect-size ranking differs across reps | One rep is an outlier -- check PCA | Drop or block as covariate; never silently include |
| No significant peaks despite obvious browser-track differences | Library-size normalization eaten the global shift | Switch to spike-in or full-library normalization |

**Operational rule:** For high-confidence reporting, require concordant detection in two methods from different families (DiffBind/DESeq2-style on consensus peaks AND csaw-style sliding windows agreeing within +/- 500 bp). Report the intersection as primary; the union as exploratory.

## Effect Size and Threshold Selection

| Question | Threshold | Rationale |
|----------|-----------|-----------|
| Statistical significance | FDR < 0.05 | Standard BH FDR (DESeq2 / edgeR / DiffBind default) |
| Stringent biological change | abs(log2FC) >= 1 (= 2-fold) | Within-noise effects below 2-fold are unreliable in chromatin |
| Conservative reporting | FDR < 0.01 AND abs(log2FC) >= 1 | Per ENCODE differential reporting guidance |
| Exploratory / discovery | FDR < 0.1 OR shrunken log2FC >= 0.585 (1.5x) | For follow-up validation, not final claim |
| Proper effect-size reporting | Use shrunken log2FC (apeglm or DESeq2 lfcShrink) when n >= 3 | Raw log2FC at low counts is volatile |

abs(log2FC) >= 1 is *not* universal. ATAC effects in primary cells (immune subsets, neurons) often max at 1.5-fold; require log2FC >= 0.585 with FDR < 0.05 for those settings.

## Hidden Batch with SVA / RUVseq

```r
library(DESeq2); library(sva)

dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
dds <- estimateSizeFactors(dds)
dat <- counts(dds, normalized=TRUE)
dat <- dat[rowMeans(dat) > 1, ]

mod  <- model.matrix(~condition, colData(dds))
mod0 <- model.matrix(~1, colData(dds))
svobj <- svaseq(dat, mod, mod0, n.sv=2)

dds$SV1 <- svobj$sv[, 1]; dds$SV2 <- svobj$sv[, 2]
design(dds) <- ~SV1 + SV2 + condition
dds <- DESeq(dds)
res <- results(dds, contrast=c('condition', 'treated', 'control'))
```

RUVseq is the alternative when negative-control regions (ChrM peaks NOT changing) or technical replicates are available. SVA is preferred when no controls exist.

## Spike-in Normalization (Reske 2020)

**Trigger:** Treatment causes whole-genome accessibility shift (HDAC inhibitor, DNMT inhibitor); RPM/CPM/RiP normalization erases the global biology.

**Mechanism:** Exogenous chromatin spike-in (Drosophila S2 nuclei) is added at constant cell number ratio pre-Tn5; reads aligning to dm6 quantify the constant exogenous baseline. Sample-level scaling factor = inverse of dm6 reads per sample, applied to human-aligned counts.

```r
library(DESeq2)

# spike_counts: per-sample dm6 read counts (one column per sample)
sf_spike <- 1 / spike_counts
sf_spike <- sf_spike / mean(sf_spike)                 # Geometric mean = 1 for stability

dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
sizeFactors(dds) <- sf_spike                          # Override library-size factors
dds <- DESeq(dds)
res <- results(dds, contrast=c('condition', 'treated', 'control'))
```

After spike-in normalization, log2FC reflects absolute accessibility change (not just relative redistribution). ENCODE 4 increasingly recommends spike-in for global-shift biology.

## Permutation Testing for Low Replicate Designs

**Trigger:** n=2 per condition; parametric NB tests give over-confident p-values.

**Mechanism:** csaw provides a permutation framework: the null is generated by shuffling sample labels; test statistic is the count-difference per window; per-region p is the rank under permutation.

```r
library(csaw); library(edgeR)

# Standard csaw counts (windows or peaks)
counts <- regionCounts(bam_files, regions, ext=200)

# Standard NB fit
y <- DGEList(counts=assay(counts), group=condition)
y <- calcNormFactors(y, method='TMM')
design <- model.matrix(~condition)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# Permutation: shuffle group labels n_perms times; track per-region rank statistic
n_perms <- 1000
perm_p <- replicate(n_perms, {
    shuffled <- sample(condition)
    design_p <- model.matrix(~shuffled)
    fit_p <- glmQLFit(estimateDisp(y, design_p), design_p)
    glmQLFTest(fit_p, coef=2)$table$F
})
observed_F <- glmQLFTest(fit, coef=2)$table$F
permp <- rowMeans(perm_p >= observed_F)
```

Permutation requires ~1000 shuffles for stable per-region p; computationally expensive but essential when parametric tests cannot be trusted.

## DESeq2 Likelihood Ratio Test for Time-Courses

```r
library(DESeq2); library(splines)

# Spline-modeled time course (5+ time points)
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata,
                              design=~ns(timepoint, df=3) + condition + ns(timepoint, df=3):condition)
dds_full <- DESeq(dds, test='LRT', reduced=~ns(timepoint, df=3) + condition)
res <- results(dds_full)
```

The LRT compares the full model (with time:condition interaction) to a reduced model without; significant peaks have time-dependent condition response. Use `df=3` natural splines for typical 5-7 timepoints; df=4-5 for >= 8.

## Hi-C-Loop-Anchored Differential

**Trigger:** Combined ATAC-seq + Hi-C/HiChIP datasets; want to test enhancer-promoter pair-level differential.

**Mechanism:** Aggregate peak-level differential signal at HiCCUPS loop anchors (or ABC-predicted enhancer-gene pairs). Combined enhancer + promoter accessibility change has more statistical power than either alone.

```r
# Pseudo-pattern: per loop, sum DESeq2 log2FC at both anchors
loops <- import('hiccups_loops.bedpe', format='bedpe')
peak_to_loop <- findOverlaps(consensus_peaks, c(anchorOne(loops), anchorTwo(loops)))

loop_lfc <- aggregate(res$log2FoldChange[queryHits(peak_to_loop)],
                      by=list(loop=ceiling(subjectHits(peak_to_loop) / 2)),
                      FUN=function(x) sum(x, na.rm=TRUE))
```

For implementation, use the `InteractionSet` Bioconductor package which preserves loop-pair structure during testing. Reference: Mumbach 2017 Nat Genet and Hwang 2024 (Hi-C loop-aggregated DA).

## Annotate Differential Peaks

```r
library(ChIPseeker); library(TxDb.Hsapiens.UCSC.hg38.knownGene)

diff_peaks <- dba.report(dba)
peakAnno <- annotatePeak(diff_peaks, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
                         tssRegion=c(-2000, 500), level='gene')
plotAnnoPie(peakAnno); plotDistToTSS(peakAnno)
genes <- as.data.frame(peakAnno)$geneId          # for GO enrichment via pathway-analysis/go-enrichment
```

`tssRegion=c(-2000, 500)` defines promoter as TSS-2kb to TSS+500bp; ChIPseeker default (`-3000, 3000`) over-counts promoter assignments. Adjust per cell type / organism.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| DiffBind very slow | Counting all peaks across all BAMs sequentially | `dba.count(..., bParallel=TRUE)` and provide `BPPARAM` |
| `unable to use the provided design matrix` (DESeq2) | Confounded design (e.g., batch perfectly aligns with condition) | Replicate in a way that breaks the confound, or drop the batch term |
| FDR list empty despite obvious differences | Library-size normalization removed global biology | Use spike-in or `DBA_NORM_LIB`; verify with browser tracks |
| Top peaks all on chrM | chrM not removed from BAM before counting | Always strip chrM upstream |
| dispersion estimate failure (DESeq2) | Too few peaks pass filter; too few reps | `filterByExpr` less aggressively; check rep count |
| `Error in if (any(out))`(csaw) | Window count below threshold | Reduce `bin.size`; check BAM is paired-end |
| ChIPseeker error on non-human TxDb | Wrong organism db loaded | Use `make_org_db` from biomartr or AnnotationDbi for non-model |
| Volcano plot symmetric about zero with no significant peaks | Hidden batch swamping signal | Run SVA/RUVseq |

## References

- Stark R & Brown G 2011 DiffBind (Bioconductor; canonical reference)
- Lun ATL & Smyth GK 2014 NAR 42:e95 (csaw windowed differential method)
- Love MI et al 2014 Genome Biol 15:550 (DESeq2)
- Robinson MD et al 2010 Bioinformatics 26:139 (edgeR)
- Chen Y et al 2016 F1000Res 5:1438 (edgeR-QL framework)
- Leek JT 2014 NAR 42:e161 (svaseq for hidden batch)
- Risso D et al 2014 Nat Biotechnol 32:896 (RUVseq)
- Reske JJ et al 2020 Epigenetics Chromatin 13:22 (ATAC normalization comparison; spike-in case)
- Corces MR et al 2018 Science 362:eaav1898 (iterative overlap, fixed-width 501 bp consensus)
- Yu G et al 2015 Bioinformatics 31:2382 (ChIPseeker)

## Related Skills

- atac-seq/atac-peak-calling - Generate per-replicate peaks
- atac-seq/consensus-peakset - Build the differential-ready consensus peakset
- atac-seq/atac-qc - Pre-screen and drop failing replicates
- atac-seq/single-cell-atac - Pseudobulk-level differential per cluster
- atac-seq/co-accessibility - Identify cis-regulatory connections among DA peaks
- differential-expression/deseq2-basics - Underlying DESeq2 patterns
- differential-expression/de-results - Effect-size reporting and shrinkage
- chip-seq/differential-binding - Same DiffBind workflow, ChIP context
- pathway-analysis/go-enrichment - Downstream gene-level enrichment of DA-associated genes
