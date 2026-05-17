---
name: bio-chipseq-differential-binding
description: Identifies differentially bound ChIP-seq regions between conditions using DiffBind, csaw (sliding windows), DESeq2/edgeR/PyDESeq2 on count matrices, NormR (control-aware), or MAnorm2. Distinguishes three distinct normalization problems (composition bias, trended bias, global shifts) and matches each to its appropriate fix including spike-in scaling. Use when comparing ChIP-seq binding between experimental conditions, choosing normalization for global vs local changes, integrating spike-in data, or reconciling DiffBind/DESeq2 disagreement.
tool_type: mixed
primary_tool: DiffBind
---

## Version Compatibility

Reference examples tested with: DiffBind 3.20+, DESeq2 1.42+, edgeR 4.0+, csaw 1.36+, PyDESeq2 0.4+, NormR 1.28+, MAnorm2 1.2+, ChIPseqSpikeInFree 1.6+.

DiffBind 3.0+ changed defaults: `summits=200` (was FALSE), `dba.normalize()` now required, blacklist filtering on by default, full library size normalization replaces reads-in-peaks. Always run `packageVersion('DiffBind')` and inspect `dba.normalize(obj, bRetrieve=TRUE)` to confirm what was applied.

# Differential ChIP-seq Binding

**"Compare protein-DNA binding between experimental conditions"** -> Identify regions where IP signal changes significantly, accounting for sequencing depth, composition bias, trended biases, and global shifts that confound naive normalization.

- R (BAM + peaks): `DiffBind::dba()` -> `dba.count()` -> `dba.normalize()` -> `dba.analyze()`
- R (count matrix): `DESeq2::DESeq()` or `edgeR::glmQLFTest()` on a peaks-by-samples matrix
- R (windows-based, global-shift-robust): `csaw::windowCounts()` -> `csaw::normFactors()` -> `edgeR::glmQLFTest()`
- R (control-aware): `normr::diffR(chip1.bam, chip2.bam, genome)` joint binomial mixture
- Python (count matrix): `pydeseq2.DeseqDataSet()`

Choice of normalization matters more than choice of test statistic (RLE vs TMM vs csaw bin-TMM on the same reference reads produce nearly identical results). Choose by which of the three normalization problems applies.

## The Three Distinct Normalization Problems

| Problem | Symptom on MA plot | Cause | Fix |
|---------|--------------------|-------|-----|
| **Composition bias** | Loess shifts off y=0 systematically | Few high-signal peaks dominate read counts; small fold changes look large or inverted | TMM on background 10 kb bins (csaw / DiffBind `background=TRUE`); NOT reads-in-peaks |
| **Trended bias (intensity-dependent)** | Loess curve sweeps from + to - across abundance | Library-prep efficiency varies with fragment abundance | Non-linear loess offsets (csaw `normalizeOffsets`); use cautiously — can over-normalize biology |
| **Global shift (treatment changes most peaks)** | Loess entirely shifted off y=0; mean log2FC ≠ 0 | Drug/perturbation changes the genome-wide level of binding (HDACi, BETi, EZH2i, target KD) | Spike-in scaling (ChIP-Rx); no algorithmic fix works |

**Why this matters:** Most published ChIP-seq differential analyses default to RLE/TMM on reads-in-peaks, which assumes "most peaks are unchanged." For HDAC inhibitors, BET inhibitors, EZH2 inhibitors, or any large dosage / target-knockdown experiment, this assumption is violated. The algorithm forces the median log2FC to zero, hides the real effect, and amplifies noise around the new "zero." The result can have the wrong sign.

**Diagnostic:** Plot MA loess on differential results. A loess curve that sweeps abundance indicates trended bias. A uniformly shifted loess indicates a global shift. Both patterns together indicate normalization is failing in two ways simultaneously.

## Algorithmic Taxonomy

| Tool | Treats | Statistical model | Strength | Fails when |
|------|--------|-------------------|----------|------------|
| **DiffBind 3.20+** | Consensus peaks summit ± 250 bp (default) | DESeq2 or edgeR backend | Mature; integrated counting/normalization; spike-in support; blacklist filter on | Default `summits=200` recenters peaks (wrong for broad marks); `DBA_NORM_LIB` is conservative but misses global shifts unless `background=TRUE` |
| **DESeq2 (direct)** | Predefined peaks | NB GLM | Familiar; transparent | RLE on reads-in-peaks fails for global shifts |
| **edgeR (direct)** | Predefined peaks | NB GLM with TMM; quasi-likelihood F-test | Cleaner small-sample inference; QL-F controls type-I error | TMM on peak counts unstable if peaks globally shifting |
| **csaw** (Lun & Smyth 2016) | Sliding windows (typically 150 bp width, 50 bp shift) | edgeR QL-F | Gold-standard for global shifts; bin-TMM composition bias; loess for trended biases | Slower; requires BAMs not count matrix; window-merge step adds complexity |
| **NormR** (Helmuth 2016) | Genomic bins | Binomial mixture (background + enriched) | Control-aware; identifies enrichment/depletion/background simultaneously | Bin-level not peak-level; older codebase; less integration with downstream tools |
| **MAnorm2** (Tu 2021) | Peaks | Hierarchical model with mean-variance trend | Designed for cross-condition with replicates | Less widely adopted; sparse maintenance |
| **SpikChIP** (Blanco 2021) | Peaks | Spike-in-aware | Multi-sample spike-in comparison | Niche; specific spike-in protocol assumed |
| **SpikeFlow** (2024) | End-to-end | Snakemake wrapper around DiffBind + spike-in | Automated; multiple normalization options | Wrapper; inherits experimental-design errors upstream |
| **ChIPComp** (Chen 2015) | Peaks | Joint Poisson with input | Control-aware | Older; less maintained |
| **ChIPseqSpikeInFree** (Jin 2020) | Peaks (post-hoc) | Distribution-shape inference | Detects global shift WITHOUT spike-in | Post-hoc heuristic only; not definitive; sanity check |

## Decision Tree: Choosing Normalization

| Scenario | Recommended | Tool / parameter |
|----------|-------------|------------------|
| Standard TF ChIP, local changes expected, balanced gain/loss | RLE on reads-in-peaks | DESeq2 default; DiffBind `DBA_NORM_RLE` |
| Histone marks, broad domains, local changes | TMM on background 10 kb bins | csaw `normFactors`; DiffBind `background=TRUE` |
| HDAC / BET / EZH2 inhibitor (global change) | Spike-in scaling | DiffBind `spikein=TRUE`; SpikeFlow; manual `sizeFactors()` from spike reads |
| Dosage titration, cell-cycle synchronization | Spike-in scaling | Same |
| ChIP target knockdown / degron | Spike-in or matched-input control subtraction | Spike-in preferred; bamCompare log2 ratio next |
| CUT&RUN/CUT&Tag standard | E. coli spike-in (carryover) | DiffBind custom scaling; see cut-and-run-tag |
| No spike-in available; suspect global shift | ChIPseqSpikeInFree | Post-hoc distribution-shape inference |
| Suspected trended (abundance-dependent) bias | Non-linear loess | csaw `normalizeOffsets` |
| Genome-wide enrichment/depletion analysis | NormR | Binomial mixture |
| Many conditions, large peak set | edgeR QL-F or DiffBind+edgeR | Better type-I control than DESeq2 Wald |

## Spike-In Scaling Factor Calculation

ChIP-Rx (Orlando 2014) uses Drosophila chromatin spike-in added at fixed concentration BEFORE IP. Egan 2016 protocol: 50,000 Drosophila S2 nuclei per 5M target cells.

**RRPM (reference-adjusted reads per million):**

```
scale_factor_i = min(N_spike_sample) / N_spike_sample_i
```

Apply to read counts pre-test, OR pass as `sizeFactors()` to DESeq2 / `normFactors()` to edgeR / `library.factors` to DiffBind:

```r
spike_in_reads <- c(120000, 145000, 110000, 95000)  # per-sample Drosophila read count
sample_names <- c('ctrl_1', 'ctrl_2', 'treat_1', 'treat_2')
scale_factors <- min(spike_in_reads) / spike_in_reads
names(scale_factors) <- sample_names

# DESeq2 with spike-in size factors
sizeFactors(dds) <- 1 / scale_factors  # DESeq2 expects inverse convention

# Or directly via DiffBind 3.x
dba_obj <- dba.normalize(dba_obj, spikein = TRUE,
                          library = DBA_LIBSIZE_FULL,
                          normalize = DBA_NORM_LIB)
```

**Rx-Input variant** (Fursova 2019): additionally scale by input spike-in to correct IP efficiency variation.

**Internal-control sanity check (Hammond Norris 2024 review):** After spike-in normalization, blacklist regions and constitutive housekeeping sites (U6 promoter, rRNA processing factors that are stable) should show no signal change. If they do, the normalization is broken — common causes:
- Spike-in scaling applied to peak counts instead of read counts
- Spike-in reads not deduplicated before scaling
- Spike-in genome not filtered for high-mapq before scaling
- Spike-in saturated (always 100k+ reads); check titration linearity

Per Hammond Norris 2024, ~25% of published spike-in ChIP papers have one of these errors.

## DiffBind Workflow (BAMs + Peaks)

**Goal:** Run DiffBind from a sample sheet to consensus peaks, normalized counts, and tested differential binding.

**Approach:** Build sample sheet, count reads in consensus peaks (summit-centered for narrow marks, full peak width for broad), choose normalization based on the three-problem framework above, then test with DESeq2 or edgeR backend.

```r
library(DiffBind)

# Sample sheet: SampleID, Condition, Replicate, bamReads, bamControl, Peaks, PeakCaller
dba_obj <- dba(sampleSheet = 'samples.csv')

# Counting: summits=250 (narrow); FALSE for broad histones (use full peak width)
dba_obj <- dba.count(dba_obj, summits = 250, minOverlap = 2, bParallel = TRUE)

# Normalization — choose per the three-problem framework
dba_obj <- dba.normalize(dba_obj,                      # default: full library size
                          method = DBA_DESEQ2,
                          normalize = DBA_NORM_LIB,
                          library = DBA_LIBSIZE_FULL)

# Background bin TMM for composition bias / broad marks
# dba_obj <- dba.normalize(dba_obj, background = TRUE)

# Spike-in scaling for global shifts
# dba_obj <- dba.normalize(dba_obj, spikein = TRUE)

dba_obj <- dba.contrast(dba_obj, design = '~ Condition')
dba_obj <- dba.analyze(dba_obj, method = DBA_DESEQ2)

# Always inspect what was actually applied
dba.normalize(dba_obj, bRetrieve = TRUE)
```

DiffBind 3.20+ defaults (verified via `bRetrieve`):
- `summits = 200` — narrow recentering (set `FALSE` for broad histones)
- `library = DBA_LIBSIZE_FULL` — full library, conservative
- `normalize = DBA_NORM_LIB` — library size only, not reads-in-peaks RLE
- Blacklist filtering on
- `background = FALSE` — set `TRUE` for composition bias

## csaw Workflow (Windows-Based)

**Goal:** Detect differential binding from sliding windows without committing to predefined peaks; robust to composition bias via background bin TMM.

**Approach:** Count reads in overlapping windows, normalize on 10 kb bins (composition bias) and/or loess (trended bias), test with edgeR QL-F, then merge significant windows into regions.

```r
library(csaw)
library(edgeR)

bam_files <- c('ctrl_1.bam', 'ctrl_2.bam', 'treat_1.bam', 'treat_2.bam')
condition <- factor(c('ctrl', 'ctrl', 'treat', 'treat'))

# Window counts: 150 bp window, 50 bp spacing (sharp marks); 1-2 kb for broad
param <- readParam(minq = 30, pe = 'both', dedup = TRUE,
                    discard = import('hg38-blacklist.v2.bed'))
windows <- windowCounts(bam_files, width = 150, ext = 200, param = param)

# Composition bias via 10 kb bins (always applied for ChIP-seq)
bg_bins <- windowCounts(bam_files, bin = TRUE, width = 10000, param = param)
windows <- normFactors(bg_bins, se.out = windows)

# Optional: trended bias via non-linear loess (use cautiously)
# windows <- normOffsets(windows, type = 'loess')

# edgeR QL-F test
y <- asDGEList(windows)
design <- model.matrix(~condition)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)
results <- glmQLFTest(fit, coef = 2)

# Merge significant windows within 1 kb into regions
merged <- mergeResults(windows, results$table, tol = 1000, merge.args = list(max.width = 5000))
```

For broad marks, increase window width to 1-2 kb and merge tolerance to 5 kb. For TFs, 150 bp window + 50 bp spacing is standard.

## DESeq2 from Count Matrix

```r
library(DESeq2)

counts <- read.delim('counts.tsv', row.names = 1, check.names = FALSE)
coldata <- data.frame(
    condition = factor(c('ctrl', 'ctrl', 'ctrl', 'treat', 'treat', 'treat')),
    row.names = colnames(counts)
)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)

# Pre-filtering for ChIP-seq is LESS aggressive than RNA-seq (peaks already enriched)
keep <- rowSums(counts(dds)) >= 1  # remove only all-zero
dds <- dds[keep, ]
dds$condition <- relevel(dds$condition, ref = 'ctrl')

dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)  # match independent-filtering optimization

# Optional: LFC shrinkage for ranking/visualization only — does NOT change padj
library(apeglm)
resLFC <- lfcShrink(dds, coef = 'condition_treat_vs_ctrl', type = 'apeglm')
```

For spike-in normalization, set `sizeFactors(dds)` from scaling factors before `DESeq(dds)`. For background-bin TMM (composition bias) use csaw to compute size factors, then transfer.

## Per-Tool Failure Modes

### DiffBind -- `summits=200` default destroys broad marks

**Trigger:** Default DiffBind 3.x call with histone broad marks (H3K27me3, H3K9me3).

**Mechanism:** `summits=200` re-centers peaks to summit ± 200 bp, throwing away most of a 10-100 kb broad domain.

**Symptom:** Differential count for broad marks much lower than expected; signal concentrated at narrow centers of broad regions.

**Fix:** `dba.count(obj, summits = FALSE, ...)` for broad marks; OR use full-width consensus peaks; OR switch to csaw with 1-2 kb windows.

### DiffBind -- `DBA_NORM_RLE` reads-in-peaks default reverses global shifts

**Trigger:** Default (legacy DiffBind < 3.0) OR explicitly setting `normalize = DBA_NORM_RLE` on a global-shift experiment.

**Mechanism:** RLE on reads-in-peaks assumes most peaks unchanged. If 80% of peaks lose signal (e.g., EZH2 inhibitor on H3K27me3), the size factors compensate by inflating the "lost" peaks' normalized values toward control levels.

**Symptom:** Differential results show fewer peaks changed than visually obvious; or signs are wrong (gain reported where loss occurred).

**Fix:** Switch to `background = TRUE` (bin TMM); for definitive analysis use `spikein = TRUE` with ChIP-Rx.

### DESeq2 -- Pre-filtering removes condition-specific peaks

**Trigger:** Applying RNA-seq-style filter `rowSums(counts) >= 10` to ChIP-seq counts.

**Mechanism:** A peak present in treatment but absent in control has near-zero control counts. Aggressive filtering removes truly differential peaks.

**Symptom:** Significantly differential gains-of-binding peaks missing from results.

**Fix:** `rowSums(counts) >= 1` only (remove all-zero rows); accept the loss of statistical power vs. recovering true differential peaks.

### csaw -- Trended bias loess over-normalizes biology

**Trigger:** Applying `normOffsets(type='loess')` when the abundance-dependent shift IS the biology.

**Mechanism:** Loess fits a smooth curve to the MA-plot trend; if treatment uniformly increases binding at low-signal peaks (which is biology), loess interprets it as a technical trend and removes it.

**Symptom:** No differential peaks detected despite obvious treatment effect.

**Fix:** Use bin-TMM only (composition bias); apply loess only after confirming the trend is technical (e.g., it appears in IgG-only samples).

### Spike-in normalization -- Scaling factor applied to wrong layer

**Trigger:** Multiplying peak counts by spike-in scaling factor.

**Mechanism:** Spike-in factors are for read-level normalization; applied to peak counts they double-correct (peak counts already reflect mapped reads).

**Symptom:** Effect sizes shifted by 2-10× from expected biology; internal-control regions show artifactual signal change.

**Fix:** Apply spike-in via `sizeFactors(dds)` (DESeq2) or `normFactors` (edgeR) BEFORE the test; or via `bamCoverage --scaleFactor` for browser tracks. Never multiply peak-level counts.

### IDR-passing peaks not used as differential input

**Trigger:** Using per-replicate MACS calls (loose `-p 1e-2`) as DiffBind peak input.

**Mechanism:** Loose ENCODE-pattern peaks include many low-confidence calls; DiffBind's `minOverlap = 2` may not filter aggressively enough.

**Symptom:** Differential results dominated by noise at marginal peaks; high false-positive rate.

**Fix:** Pre-filter peak input to IDR-passing peaks (TF) or naive-overlap-passing peaks (histone) before DiffBind ingestion.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| DiffBind + DESeq2 directly differ wildly | Different normalization (DiffBind default = library size; DESeq2 default = RLE) | Force same normalization; differences should shrink to <5% |
| csaw windows + DiffBind peaks disagree | csaw catches sub-peak local maxima or wider regions DiffBind missed | Inspect IGV; csaw windows-based often more sensitive to broad/diffuse changes |
| Spike-in scaled + non-scaled give opposite signs | Global shift present | Spike-in is correct; non-scaled is fooled by composition bias |
| DiffBind run twice gives different results | Different `summits` or `minOverlap` settings | Verify via `dba.normalize(obj, bRetrieve=TRUE)`; pin parameters in script |
| Few replicates, large fold changes, low padj | DESeq2 dispersion estimate unstable with n=2 | Switch to edgeR QL-F or DiffBind with edgeR backend |
| Different fold changes in DiffBind 3.x vs 2.x | Default normalization changed | Match settings explicitly; document version in methods |

**Operational rule for publication-grade:** Run on ENCODE-pattern peaks (IDR or naive overlap-passing), normalize with both reads-in-peaks AND background-bin AND spike-in if available; require concordance across at least two methods. For global-shift experiments, spike-in is mandatory.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `Error: dba.normalize() must be called before dba.analyze()` | DiffBind 3.x change; not in older docs | Add `dba.normalize(obj)` step |
| DiffBind very slow on many samples | Sequential counting | `dba.count(obj, bParallel = TRUE)` |
| `Error in DESeq()` with few replicates | Dispersion estimation unstable | Use `DESeq(dds, fitType = 'parametric', sfType = 'poscounts')` or switch to edgeR |
| Spike-in size factors mostly NA | Spike-in reads not in sample sheet or wrong BAM path | Verify spike-in BAM exists; load via DiffBind `spikein` field |
| All padj = NA | Independent filtering too aggressive | Lower `alpha` in `results()` to match intended threshold |
| Volcano plot inverted (down peaks on right) | Contrast direction reversed | `contrast = c('condition', 'treat', 'ctrl')` for positive log2FC = up in treat |

## References

- Stark R & Brown G 2011 Bioconductor (DiffBind)
- Lun ATL & Smyth GK 2016 Nucleic Acids Res 44:e45 (csaw)
- Love MI et al 2014 Genome Biol 15:550 (DESeq2)
- Robinson MD et al 2010 Bioinformatics 26:139 (edgeR; TMM)
- Helmuth J et al 2016 bioRxiv (NormR)
- Tu S et al 2021 Front Genet 12:646533 (MAnorm2)
- Orlando DA et al 2014 Cell Rep 9:1163 (ChIP-Rx framework)
- Egan B et al 2016 PLoS One 11:e0166438 (ChIP-Rx protocol)
- Fursova NA et al 2019 Mol Cell 74:1020 (Rx-Input scaling)
- Jin H et al 2020 Bioinformatics 36:1270 (ChIPseqSpikeInFree)
- Blanco E et al 2021 NAR Genom Bioinform 3:lqab064 (SpikChIP)
- Hammond Norris and Norris 2024 (review of spike-in normalization failure modes; published as a PMC-indexed open-access methodology review at PMC12266361)
- Gregoricchio S et al 2024 NAR Genom Bioinform 6:lqae118 (SpikeFlow Snakemake pipeline)

## Related Skills

- chip-seq/peak-calling - Upstream peak calling for DiffBind input
- chip-seq/chipseq-qc - Replicate concordance required before differential
- chip-seq/spike-in-normalization - Detailed spike-in workflow and scaling factor calculation
- chip-seq/cut-and-run-tag - CUT&RUN/CUT&Tag uses E. coli spike-in (different from Drosophila ChIP-Rx)
- chip-seq/peak-annotation - Annotate differential peaks to genes/cCREs
- differential-expression/deseq2-basics - DESeq2 fundamentals
- differential-expression/edger-basics - edgeR quasi-likelihood framework
- atac-seq/differential-accessibility - Parallel ATAC differential workflow
