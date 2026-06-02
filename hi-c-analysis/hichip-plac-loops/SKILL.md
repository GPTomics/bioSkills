---
name: bio-hi-c-analysis-hichip-plac-loops
description: Calls significant loops from protein-directed and targeted 3C assays (HiChIP, PLAC-seq, Capture Hi-C/PCHi-C, ChIA-PET) where the contact background is peak-anchored and coverage-biased, so generic Hi-C loop callers (cooltools dots, Juicer HiCCUPS) use the wrong null. Covers FitHiChIP (config-driven coverage+distance-decay spline regression, peak-to-peak vs peak-to-all foreground, loose vs stringent background, coverage vs ICE bias), MAPS (positive Poisson regression on bias factors for PLAC-seq/HiChIP), hichipper (restriction-site-distance bias model + library QC), CHiCAGO (Delaporte two-component Brownian+technical background for asymmetric bait x other-end Capture Hi-C), the with/without separate-ChIP anchor decision, and differential loops via diffloop. Use when calling loops from HiChIP/PLAC-seq/Capture Hi-C, choosing FitHiChIP/MAPS/CHiCAGO, picking peak-to-all vs peak-to-peak, setting the loop FDR, supplying ChIP peaks as anchors, QCing a HiChIP library, or comparing loops between conditions.
tool_type: mixed
primary_tool: fithichip
---

## Version Compatibility

Reference examples tested with: FitHiChIP 11.0+, MAPS 1.1+, hichipper 0.7+, CHiCAGO 1.20+ (Bioconductor 3.18+), HiC-Pro 3.1+, diffloop 1.18+ (Bioconductor 3.18+).

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- FitHiChIP: read the shipped `configfile` comments; parameter names and defaults change between releases
- R: `packageVersion('diffloop')`/`packageVersion('Chicago')` then `?function` to check signatures

FitHiChIP is driven entirely by a key=value config file passed with `-C`; the loop caller, background, and bias model are set there, not on the command line. MAPS, hichipper, and CHiCAGO each expect a specific upstream format (HiC-Pro valid pairs / .allValidPairs, or HiCUP+capture design files for CHiCAGO). If a tool errors, introspect the installed version's config/help and adapt the example rather than retrying.

# Protein-Directed and Targeted 3C Loop Calling

**"My HiChIP/PLAC-seq/Capture Hi-C has loops anchored at CTCF/H3K27ac/promoters - which contacts are real?"** -> Call loops with a method whose null jointly models the per-anchor coverage bias AND the distance-decay, not the uniform/donut background that generic Hi-C callers assume.
- CLI: `bash FitHiChIP_HiCPro.sh -C config_fithichip` (HiChIP/PLAC-seq); edit and run `run_pipeline.sh` (MAPS; it invokes MAPS.py internally) for PLAC-seq/HiChIP
- R: `runChicago(...)` then PIRs at CHiCAGO score >= 5 (Capture Hi-C/PCHi-C); `quickAssoc()`/`loopAssoc()` (diffloop, differential)

## The Single Most Important Modern Insight -- Generic Hi-C Loop Callers Use the Wrong Null on Peak-Anchored Data

Running `cooltools dots` or Juicer HiCCUPS on HiChIP, PLAC-seq, or Capture Hi-C is a documented error, not a shortcut. Those callers test each pixel against a *local, roughly-uniform* expected background (a donut/expected neighborhood) built for a genome-wide-uniform in-situ Hi-C map. Protein-directed and capture assays violate that assumption in the most consequential way possible: the antibody (or oligo capture) **enriches contacts at the factor's binding sites**, so 1D coverage is wildly non-uniform - an H3K27ac anchor can carry 100x the read depth of a flanking non-peak bin. A donut null reads that coverage spike as contact enrichment and calls a "loop" at every peak. The dedicated callers exist precisely to fix this, and they all share one move: **regress out the per-anchor coverage bias before testing the distance-decayed contact frequency.**

Three load-bearing consequences:

1. **The hard part is 3C statistics, and it lives here; the peak-calling half lives in chip-seq.** FitHiChIP/MAPS/CHiCAGO each fit a *significance model* (spline regression on coverage + genomic distance; positive Poisson regression on bias factors; a two-component Brownian+technical background). That model - not the antibody - is the deliverable. Anchor/peak calling (where the protein binds) is chip-seq's job (-> chip-seq/peak-calling); this skill consumes those peaks and produces FDR-controlled loops.

2. **Protein-targeting buys depth efficiency, so loops are called at far lower total depth than Hi-C.** HiChIP/PLAC-seq concentrate reads onto a small anchored sub-space, needing ~5-10 read pairs per interaction versus ~100-1000 for genome-wide Hi-C (Mumbach 2016: >10x more conformation-informative reads, >100x less input than ChIA-PET). A 100-200M-pair HiChIP library calls loops that would need billions of pairs in Hi-C - but only at the protein's anchors, and only with the right null.

3. **"Peaks from the same data" is a circularity trap.** When no separate ChIP-seq exists, HiChIP-derived peaks (hichipper, HiChIP-Peaks) are used as anchors - but calling peaks and loops from the same reads couples the two error structures. Prefer an independent ChIP-seq peak set as the anchor reference when one exists; if not, use a HiChIP-native peak caller and treat anchor confidence as part of the loop's uncertainty, not a given.

## Method Taxonomy

| Tool | Assay | Null / significance model | Anchors | When |
|------|-------|---------------------------|---------|------|
| FitHiChIP | HiChIP, PLAC-seq, (CHi-C, ChIA-PET) | spline regression of contact count on coverage bias AND genomic distance; loose (peak-to-all) vs stringent (peak-to-peak) background; coverage-bias or ICE-bias regression | ChIP/HiChIP peak file | default; recovers Hi-C/CHi-C/ChIA-PET contacts best (Bhattacharyya 2019); config-driven |
| MAPS | PLAC-seq, HiChIP | zero-truncated (positive) Poisson regression removing effective-fragment/GC/mappability AND ChIP-enrichment bias, then test normalized frequency at anchored bins | AND-set vs XOR-set anchored bins | model-based PLAC-seq/HiChIP, 4DN-adopted; two-step (bias model -> significance) |
| hichipper | HiChIP | background read density modeled as a function of proximity to restriction sites; loop strength + confidence per anchor | self-derived (restriction-aware) | restriction-aware QC + loop calling without separate ChIP; feeds diffloop |
| CHiCAGO | Capture Hi-C, PCHi-C | Delaporte two-component background: Brownian (distance-dependent, NB) + technical (distance-independent, Poisson), fit per bait; report PIRs at score >= 5 | baited fragments (asymmetric bait x other-end) | promoter/region-capture; asymmetric design where dots cannot apply |
| HiChIP-Peaks | HiChIP | peak calling from HiChIP signal (not loop calling) | n/a (produces anchors) | when no separate ChIP exists and anchors must come from HiChIP itself |
| diffloop | any loop set (HiChIP/ChIA-PET) | edgeR-style count test on a union loop set across conditions | from the union set | differential looping between conditions (not a caller) |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| HiChIP/PLAC-seq, have a separate ChIP-seq peak set | FitHiChIP with `PeakFile=` the ChIP peaks | independent anchors break the peak/loop circularity; FitHiChIP is the default |
| PLAC-seq/HiChIP, prefer a regression-model caller | MAPS | positive Poisson regression explicitly removes ChIP-enrichment bias |
| No separate ChIP; need anchors + library QC fast | hichipper (then FitHiChIP/diffloop) | restriction-aware, self-derives anchors, reports library quality |
| Capture Hi-C / Promoter-Capture Hi-C | CHiCAGO, PIRs at score >= 5 | asymmetric bait x other-end; two-component per-bait background |
| H3K27ac/broad anchors, want sensitivity | FitHiChIP loose (peak-to-all) background | most contacts have at least one peak anchor |
| CTCF/cohesin sharp anchors, want specificity | FitHiChIP stringent (peak-to-peak) background | restricts foreground to peak-peak contacts |
| Compare loops between conditions | -> hic-differential context; quantify with diffloop / FitHiChIP DiffAnalysis | union anchors + count test, not pixel subtraction |
| Generic in-situ Hi-C (no protein/capture) | -> loop-calling (cooltools dots / chromosight) | uniform background is correct there; do NOT use it here |
| Anchors not yet called | -> chip-seq/peak-calling | peak calling is chip-seq's competency; this skill consumes peaks |
| Annotate loop anchors with TFs/genes | -> chip-seq/peak-annotation, atac-seq/enhancer-gene-linking | anchor-to-feature assignment lives there |

## FitHiChIP - Coverage + Distance Spline Regression (default)

**Goal:** Call FDR-controlled loops from a HiChIP/PLAC-seq library whose anchors are defined by an (ideally independent) ChIP-seq peak set.

**Approach:** FitHiChIP reads valid pairs (HiC-Pro format), bins them, fits a spline regression of contact count on BOTH the genomic-distance decay and the per-bin coverage bias, then assigns each candidate contact an FDR. Everything - resolution, foreground type, background, bias model, FDR - is set in a key=value config passed with `-C`; the command line itself takes no analysis parameters.

```bash
# config_fithichip (key=value; comments stripped). Run: bash FitHiChIP_HiCPro.sh -C config_fithichip
ValidPairs=sample.allValidPairs.gz   # HiC-Pro valid pairs (or set Matrix=/Interval= for matrix input)
PeakFile=chipseq_peaks.bed           # anchors; prefer an INDEPENDENT ChIP-seq peak set over HiChIP-derived
ChrSizeFile=hg38.chrom.sizes
OutDir=fithichip_out/
PREFIX=sample
BINSIZE=5000                          # 5kb: standard HiChIP anchor resolution (~2.5kb effective, hichipper)
LowDistThr=20000                      # 20kb floor: below this, contacts are dominated by self-ligation/diagonal
UppDistThr=2000000                    # 2Mb ceiling: loops beyond this are rare and noise-dominated
IntType=3                             # 3=peak-to-all (loose foreground); 1=peak-to-peak (stringent)
UseP2PBackgrnd=0                      # 0=loose (peak-to-all) background; 1=stringent (peak-to-peak) background
BiasType=1                            # 1=coverage-bias regression (default); 2=ICE-bias regression
MergeInt=1                            # merge adjacent significant contacts into one loop (recommended)
QVALUE=0.01                           # FDR cutoff for significant loops
```

`IntType` sets the foreground (which candidate contacts are tested); `UseP2PBackgrnd` sets the background the regression is fit against. The two together encode the loose-vs-stringent choice: peak-to-all foreground + loose background maximizes sensitivity for broad marks (H3K27ac); peak-to-peak foreground + stringent background maximizes specificity for sharp factors (CTCF/cohesin). Output significant loops land under a nested `OutDir/FitHiChIP_Peak2ALL_b<bin>_L<low>_U<upp>/P2Pbckgr_<0|1>/.../` tree, in `<PREFIX>.interactions_FitHiC_Q<QVALUE>.bed` (and `..._MergeNearContacts.bed` when `MergeInt=1`); locate it with `find OutDir -name '*interactions_FitHiC_Q*.bed'`.

## MAPS - Positive Poisson Regression on Bias Factors

**Goal:** Call PLAC-seq/HiChIP loops with an explicit regression model that removes both the generic 3C biases and the ChIP-enrichment bias.

**Approach:** MAPS is a two-step pipeline: first fit a zero-truncated (positive) Poisson regression of observed contact counts on effective-fragment length, GC content, mappability, and ChIP-enrichment per bin; then test each anchored bin-pair's count against the model-normalized expectation, controlling FDR. It distinguishes AND anchors (both ends in a peak) from XOR anchors (one end), reflecting the peak-to-peak vs peak-to-all distinction.

```bash
# MAPS is driven by a COPIED run_pipeline.sh with key=value bash variables, not CLI flags.
# Edit run_pipeline_sample.sh, then run it: ./run_pipeline_sample.sh
bin_size=5000                              # 5kb anchor bin
binning_range=1000000                      # max interaction distance modeled
fdr=2                                       # -log10(FDR) cutoff; 2 means FDR <= 0.01
dataset_name='sample'
macs2_filepath='chipseq_peaks.narrowPeak'  # ChIP/HiChIP anchors
organism='hg38'                            # selects the bundled effective-length/GC/mappability bias track
# run_pipeline.sh runs feather (preprocessing) then MAPS.py (positive Poisson regression) internally
```

The model-based design is MAPS's signature: it does not subtract a local background; it predicts each bin-pair's expected count from the bias covariates and flags positive residuals. FitHiChIP and MAPS disagree substantially on the same data (the literature reports tens-of-thousands-loop differences at matched FDR) - the model assumptions differ, so report which caller and its settings.

## CHiCAGO - Two-Component Background for Capture Hi-C

**Goal:** Call significant promoter-interacting regions (PIRs) from Capture Hi-C / PCHi-C, where oligo capture makes the map asymmetric (baited fragment x any other-end).

**Approach:** CHiCAGO fits a per-bait background with two components - a Brownian (distance-dependent, negative-binomial) term and a technical-noise (distance-independent, Poisson) term, convolved as a Delaporte distribution - then scores each bait-other-end pair as a weighted -log p-value; report other-ends above the conventional score threshold.

```r
# Reference: Chicago 1.20+ (Bioconductor 3.18+) | Verify API if version differs
library(Chicago)

CHICAGO_SCORE <- 5   # conventional PIR threshold (Cairns 2016); soft 3-5 grey zone, >=5 = called
cd <- setExperiment(designDir = 'capture_design/')         # baitmap/rmap/NPB/NBaitsPB/proxOE from the capture design
cd <- readAndMerge(files = c('sample_rep1.chinput', 'sample_rep2.chinput'), cd = cd)
cd <- chicagoPipeline(cd)                                  # fits Brownian+technical background, scores all bait x other-end
exportResults(cd, file.path('chicago_out', 'sample'), format = 'washU_text')   # PIRs at score >= CHICAGO_SCORE
```

Neither cooltools dots nor FitHiChIP's symmetric model applies to Capture Hi-C: the bait-vs-other-end asymmetry and the per-bait normalization are the whole point. The capture design files (`baitmap`, `rmap`, and the precomputed `NPB`/`NBaitsPB`/`proxOE` from `makeDesignFiles.py`) encode which fragments were baited and the distance-binned background normalization.

## Differential Loops with diffloop

**Goal:** Find loops that change strength between conditions, given per-condition loop call sets.

**Approach:** Build a UNION loop set across all samples, count the read pairs supporting each loop per replicate, then run an edgeR-style count test on the union set; there is no "DESeq2 for loops," so the union-then-count workflow is the standard.

```r
# Reference: diffloop 1.18+ (Bioconductor 3.18+) | Verify API if version differs
library(diffloop)

loops <- loopsMake(beddir = 'hichipper_loops/')            # reads the hichipper-preprocessed loop directory
loops <- subsetLoops(loops, loops@rowData$loopWidth >= 20000)   # drop sub-20kb (self-ligation regime)
groups <- c('wt', 'wt', 'ko', 'ko')
loops <- updateLDGroups(loops, groups)
res <- quickAssoc(loops)   # two-group edgeR exact test on the union set; loopAssoc(loops, coef=, design=) for a GLM
```

diffloop pairs naturally with hichipper output. FitHiChIP also ships a differential-analysis script (`DiffAnalysisHiChIP.r`); either way the unit of comparison is a union anchor/loop set, not a per-pixel matrix subtraction (-> hic-differential for the matrix-level framing).

## Per-Method Failure Modes

### Generic dots/HiCCUPS on protein-directed data
**Trigger:** running `cooltools dots` or Juicer HiCCUPS on a HiChIP/PLAC-seq/Capture cooler. **Mechanism:** the donut/local-expected null assumes uniform coverage; antibody/capture enrichment spikes coverage at anchors. **Symptom:** a "loop" at essentially every peak; calls that do not reproduce across replicates. **Fix:** use FitHiChIP/MAPS (HiChIP/PLAC-seq) or CHiCAGO (Capture); they regress out coverage bias.

### Anchors and loops from the same reads (circularity)
**Trigger:** HiChIP-derived peaks used as the FitHiChIP `PeakFile` when a separate ChIP-seq exists. **Mechanism:** peak and loop errors share a source, inflating apparent confidence at high-coverage anchors. **Symptom:** loops concentrate at the strongest coverage peaks regardless of biology. **Fix:** anchor on an independent ChIP-seq peak set; reserve HiChIP-native peaks for when no ChIP exists.

### Wrong foreground/background pair for the mark
**Trigger:** stringent peak-to-peak background on a broad H3K27ac library, or loose on sharp CTCF. **Mechanism:** the foreground/background must match the anchor sharpness. **Symptom:** too few loops (over-stringent on broad marks) or noisy excess (over-loose on sharp factors). **Fix:** loose/peak-to-all for broad marks; stringent/peak-to-peak for CTCF/cohesin.

### CHiCAGO design files mismatched to the capture
**Trigger:** baitmap/rmap or precomputed NPB/NBaitsPB/proxOE built for a different fragmentation or bait set. **Mechanism:** the per-bait background normalization depends on the exact capture design. **Symptom:** absurd scores or empty PIR lists. **Fix:** regenerate design files with `makeDesignFiles.py` from the actual rmap/baitmap.

### Short-range contacts not excluded
**Trigger:** `LowDistThr` left at 0 / no distance floor. **Mechanism:** sub-20kb contacts are dominated by the diagonal, self-ligation, and re-ligation, not loops. **Symptom:** a wall of "loops" hugging the diagonal. **Fix:** set a distance floor (FitHiChIP `LowDistThr=20000`; equivalent in MAPS/CHiCAGO).

### Chromosome-name mismatch across inputs
**Trigger:** `chr1` in the valid pairs vs `1` in the peak/chrom-size file. **Mechanism:** anchors silently fail to intersect the contacts. **Symptom:** few or zero loops, no error. **Fix:** harmonize chromosome naming across valid pairs, PeakFile, and ChrSizeFile.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Loop resolution 5kb (HiChIP) | hichipper effective ~2.5kb (Lareau & Aryee 2018) | standard HiChIP anchor bin; finer needs more depth, coarser blurs anchors |
| Lower distance floor 20kb | FitHiChIP default; self-ligation/diagonal regime | below ~20kb contacts are dominated by religation/dangling/diagonal, not loops |
| Upper distance ceiling 2Mb | FitHiChIP default | loops beyond ~2Mb are rare and noise-dominated at typical HiChIP depth |
| Loop FDR (q) <= 0.01 | FitHiChIP/MAPS default | genome-wide candidate-contact testing needs strict FDR; 0.05 acceptable for discovery |
| CHiCAGO PIR score >= 5 | Cairns 2016 convention | weighted -log p threshold; 3-5 is a soft grey zone, >=5 is called |
| ~5-10 read pairs per interaction | Mumbach 2016 (HiChIP efficiency) | protein-targeting lets loops be called at far lower depth than Hi-C |
| MAPS/FitHiChIP at 5-10kb, ~100-300M valid pairs | HiChIP depth practice | anchored sub-space is small, so usable loop resolution arrives well below Hi-C billions |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| A loop at every peak; no replicate reproducibility | generic dots/HiCCUPS used on HiChIP/capture | switch to FitHiChIP/MAPS/CHiCAGO |
| FitHiChIP runs but finds almost nothing | over-stringent background on a broad mark, or wrong `PeakFile` | use loose/peak-to-all; verify the peak set matches the antibody |
| FitHiChIP `PeakFile`/`ChrSizeFile` error | missing mandatory config key or wrong path | every mandatory key (`PeakFile`, `ChrSizeFile`, `OutDir`) must be set |
| Few/zero loops, no error | chrom-name mismatch (`chr1` vs `1`) across inputs | harmonize naming across valid pairs, peaks, chrom sizes |
| CHiCAGO empty/absurd PIR list | design files mismatched to the capture | regenerate baitmap/rmap + NPB/NBaitsPB/proxOE with `makeDesignFiles.py` |
| diffloop `loopsMake` reads nothing | wrong bedpe directory or per-sample naming | point `beddir` at the per-sample hichipper loop bedpe files |
| Wall of diagonal-hugging loops | no lower distance threshold | set `LowDistThr` (FitHiChIP) / equivalent distance floor |

## References

- Mumbach MR, Rubin AJ, Flynn RA, Dai C, Khavari PA, Greenleaf WJ, Chang HY. 2016. HiChIP: efficient and sensitive analysis of protein-directed genome architecture. *Nat Methods* 13(11):919-922.
- Fang R, Yu M, Li G, Chee S, Liu T, Schmitt AD, Ren B. 2016. Mapping of long-range chromatin interactions by proximity ligation-assisted ChIP-seq. *Cell Res* 26:1345-1348.
- Bhattacharyya S, Chandra V, Vijayanand P, Ay F. 2019. Identification of significant chromatin contacts from HiChIP data by FitHiChIP. *Nat Commun* 10:4221.
- Juric I, Yu M, Abnousi A, Raviram R, Fang R, Zhao Y, Zhang Y, Qiu Y, Hu M, et al. 2019. MAPS: model-based analysis of long-range chromatin interactions from PLAC-seq and HiChIP experiments. *PLoS Comput Biol* 15(4):e1006982.
- Lareau CA, Aryee MJ. 2018. hichipper: a preprocessing pipeline for calling DNA loops from HiChIP data. *Nat Methods* 15:155-156.
- Cairns J, Freire-Pritchett P, Wingett SW, Varnai C, Dimond A, Plagnol V, Zerbino D, Schoenfelder S, Javierre BM, Osborne C, Fraser P, Spivakov M. 2016. CHiCAGO: robust detection of DNA looping interactions in Capture Hi-C data. *Genome Biol* 17:127.
- Lareau CA, Aryee MJ. 2018. diffloop: a computational framework for identifying and analyzing differential DNA loops from sequencing data. *Bioinformatics* 34(4):672-674.

## Related Skills

- loop-calling - The bulk in-situ Hi-C counterpart (cooltools dots / chromosight); correct null there, wrong null here
- hic-differential - Matrix-level condition comparison framing behind differential loops
- contact-pairs - Produces the valid pairs (HiC-Pro / pairtools) these callers consume
- hic-data-io - Cooler handling of the contact maps upstream of anchored loop calling
- chip-seq/peak-calling - Calls the ChIP-seq peaks used as independent loop anchors
- chip-seq/peak-annotation - Annotate loop anchors with TFs/genes
- atac-seq/enhancer-gene-linking - Enhancer-promoter contacts complement HiChIP/PCHi-C loops
- genome-intervals/overlap-significance - Test loop-anchor enrichment at features against a structured null
