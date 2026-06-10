---
name: bio-methylation-dmr-detection
description: Detects differentially methylated regions (DMRs) from short-read bisulfite (WGBS/RRBS), array, and long-read methylation count tables using dmrseq (permutation region-FDR over the region selection), DSS callDMR (beta-binomial), methylKit tiles, bsseq BSmooth, DMRcate Gaussian-kernel smoothing, metilene, and comb-p. Covers why a DMR is DEFINED by arbitrary thresholds (min-CpGs, max-gap, delta-beta, q) and a smoothing bandwidth, why selecting extreme runs of CpGs then testing them on the same data is post-selection inference, why region q-values are not comparable across tools, and a single-sample domain-segmentation section (PMD, UMR/LMR, MethylSeekR, solo-WCGW) that must run before focal calling on cancer/aging genomes. Use when calling region-level methylation differences, choosing a DMR caller, controlling region-level FDR, or segmenting megabase methylation domains. For per-site testing see differential-cpg-testing; for the methylKit object model see methylkit-analysis.
tool_type: r
primary_tool: dmrseq
---

## Version Compatibility

Reference examples tested with: dmrseq 1.22+, DSS 2.50+, methylKit 1.28+, bsseq 1.38+, DMRcate 2.16+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The GENOME BUILD is a version that matters. methylKit/bsseq/DMRcate `assembly=` is metadata, but annotation packages (`annotatr::build_annotations(genome='hg38')`, `TxDb.Hsapiens.UCSC.hg38.knownGene`) are build-specific and must match the alignment genome. DMRcate `arraytype='EPIC'`/`'450K'` and the IlluminaHumanMethylation annotation package set the array CpG universe. DMRcate defaults shift across Bioconductor releases (the `C` kernel scaling has no single default - it is platform-dependent); confirm with `?dmrcate` on the installed build.

# DMR Detection

**"Find differentially methylated regions"** -> Detect candidate regions, then score them with a null that re-ran the selection - because a DMR is defined by a chain of thresholds, not found, and only a selection-aware q-value is honest.
- R: `dmrseq(bs, testCovariate='condition')` (selection-aware) ; `DSS::callDMR()`, `methylKit::tileMethylCounts()`, `bsseq` BSmooth, `DMRcate::dmrcate()` (combine-and-correct)

Scope: REGION-level differential methylation from any per-CpG methylation+coverage table (WGBS/RRBS counts, array beta/M, long-read modkit bedMethyl), plus single-sample domain segmentation. Per-site DMC/DMP testing -> differential-cpg-testing. The methylKit import/filter/unite object model -> methylkit-analysis. Long-read MM/ML calling that produces the counts -> long-read-sequencing/nanopore-methylation. Functional enrichment of DMR genes -> pathway-analysis/go-enrichment (with the CpG-bias correction noted below).

## The Single Most Important Modern Insight -- A DMR Is DEFINED, Not Found, and the Region p-Value Is Only Honest If Its Null Re-Ran the Selection

The naive recipe - compute a per-CpG statistic, select runs of CpGs that look extreme to DEFINE candidate regions, then test those same regions and report a p/q on the same data - reuses the data twice. The regions were CHOSEN because they were extreme; scoring them with the data that selected them inflates significance and produces uncalibrated region FDR. This is post-selection inference, the field's original sin. No error is thrown; the q-values just lie. Three corollaries:

1. **dmrseq neutralizes the selection; the others combine-and-correct but do not re-select.** dmrseq (Korthauer 2019 *Biostatistics* 20:367) builds the null by PERMUTING the condition labels and RE-RUNNING the entire candidate-detection procedure on each permutation, pooling permuted region statistics into one genome-wide null - so its q-value is on REGIONS and accounts for selection. comb-p and DMRcate combine existing per-CpG p-values (Stouffer/Fisher/SLK) and apply a region multiplicity correction, but never re-run selection. methylKit tiles use FIXED windows (boundaries not chosen from the data, so the selection part is sidestepped) but ignore inter-tile correlation. DSS callDMR merges significant CpGs (threshold-then-define).

2. **Region q-values are NOT comparable across tools.** "I found 4,000 DMRs at q<0.05" is meaningless without naming the tool and what its q controls. Cross-tool OVERLAP is the real evidence - run two callers and intersect.

3. **Thresholds are conventions, not biology.** delta-beta 25%, min-CpGs 3, max-gap 1000bp, q<0.01 are tutorial folklore. The same data yields radically different DMR sets under different settings. Report and justify every knob; never present one config as correct.

Organize the analysis around defending these, not around listing functions.

## Tool Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| dmrseq | Korthauer 2019 *Biostatistics* 20:367 | GLS area statistic + permutation null that re-runs selection; smooths the difference internally | headline WGBS inference; calibrated region FDR; >=2 reps/group |
| DSS callDMR | Feng 2014 *Nucleic Acids Res* 42:e69; Park & Wu 2016 *Bioinformatics* 32:1446 | per-CpG Bayesian beta-binomial dispersion shrinkage, then merge significant CpGs | small n, complex/multi-factor designs; low-coverage with smoothing |
| methylKit tiles | Akalin 2012 *Genome Biol* 13:R87 | fixed windows + logistic/F test per tile | fast RRBS/WGBS screening; reuses the methylKit object model |
| bsseq (BSmooth) | Hansen 2012 *Genome Biol* 13:R83 | per-sample local-likelihood smoothing -> smoothed t-statistic | low/uneven-coverage WGBS; superseded by dmrseq for calibrated FDR |
| DMRcate | Peters 2015 *Epigenetics Chromatin* 8:6; Peters 2021 *Nucleic Acids Res* 49:e109 | Gaussian-kernel smoothing of the per-CpG statistic | arrays (450K/EPIC) and WGBS (different kernel) |
| metilene | Juhling 2016 *Genome Res* 26:256 | binary segmentation + 2D Kolmogorov-Smirnov; standalone C CLI | fast whole-genome second caller; data-adaptive boundaries; tolerates missingness |
| comb-p | Pedersen 2012 *Bioinformatics* 28:2986 | ACF + Stouffer-Liptak-Kechris combination + Sidak; Python CLI on any p-values | array EWAS regions; tool-agnostic corroboration of any per-CpG p-values |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| WGBS, >=2 reps/group, headline region inference | dmrseq | only caller whose region FDR accounts for selection |
| WGBS, small n, multi-factor / covariates | DSS (`DMLtest.multiFactor` -> `callDMR`) | beta-binomial dispersion shrinkage + model formula |
| Low / uneven coverage WGBS | dmrseq (smooths internally) or DSS `smoothing=TRUE` | smoothing borrows strength across CpGs |
| RRBS, quick region screen | methylKit tiles (filtered, `cov.bases>=3`) | fixed windows; fast; CpG-island enriched; screening q, not selection-corrected |
| EPIC/450K array | DMRcate array mode (`cpg.annotate('array')`) or comb-p on limma p | kernel tuned for array spacing; methylKit/bsseq/dmrseq are count-based |
| Corroborate any DMR set | run a second caller, intersect | region q is not comparable across tools |
| Cancer / aging / placenta / cultured-cell WGBS | segment PMDs FIRST (see domain section) | a focal caller manufactures fake hypo-DMRs from PMD background |
| Per-CpG, not region | -> differential-cpg-testing | site-level test before aggregating to regions |
| Long-read modkit bedMethyl input | -> long-read-sequencing/nanopore-methylation, then any caller here | modkit calls; feed Nmod + Nvalid_cov here for region statistics |

## dmrseq: The Selection-Aware Headline Caller

**Goal:** Call WGBS DMRs with a region-level FDR that survives the region-selection step.

**Approach:** Build a `BSseq` object from raw counts (do NOT pre-smooth), filter loci with zero coverage in any sample, then run `dmrseq`, which detects candidate regions (runs of CpGs whose smoothed methylation-difference coefficient exceeds `cutoff`), computes one GLS area statistic per region, and generates the null by permuting labels and re-running detection.

```r
library(dmrseq)
library(bsseq)

bs <- read.bismark(c('ctrl1.cov.gz', 'ctrl2.cov.gz', 'treat1.cov.gz', 'treat2.cov.gz'),
                   colData = DataFrame(condition = c('ctrl', 'ctrl', 'treat', 'treat')),
                   rmZeroCov = TRUE, strandCollapse = TRUE)

# dmrseq requires every locus to have non-zero coverage in EVERY sample.
bs <- bs[rowSums(getCoverage(bs) == 0) == 0, ]

# cutoff=0.1 only SEEDS candidate detection (10% smoothed difference); significance
# comes from the permutation statistic, so dmrseq does NOT hard-threshold delta-beta.
dmrs <- dmrseq(bs, testCovariate = 'condition', cutoff = 0.1)
sig <- dmrs[dmrs$qval < 0.05]   # qval is a REGION FDR that accounts for selection
```

dmrseq smooths the DIFFERENCE internally (`bpSpan`/`minInSpan`/`maxGapSmooth`); running `BSmooth()` first double-smooths and invalidates the model. The permutation null is COARSE at 2-vs-2 (few distinct label permutations) - the package pools across all candidate regions to compensate, but more replicates give a finer null. Use `adjustCovariate` for nuisance variables and `block = TRUE` for large-scale differential blocks.

## DSS: Beta-Binomial Dispersion Shrinkage

**Goal:** Call DMRs with per-CpG dispersion shrinkage for small n or a multi-factor design.

**Approach:** `DMLtest` does per-CpG Wald tests with a Bayesian beta-binomial dispersion estimate, then `callDMR` merges significant CpGs into regions.

```r
library(DSS)

bs <- makeBSseqData(list(c1, c2, t1, t2), c('C1', 'C2', 'T1', 'T2'))   # each: chr/pos/N/X data.frame
dml <- DMLtest(bs, group1 = c('C1', 'C2'), group2 = c('T1', 'T2'), smoothing = TRUE)   # TRUE for low-cov WGBS

# callDMR defaults (verify on installed build): delta=0, p.threshold=1e-5,
# minlen=50, minCG=3, dis.merge=100, pct.sig=0.5.
# delta=0 means NO effect-size floor - SET it explicitly so tiny shifts are not called.
dmrs <- callDMR(dml, delta = 0.1, p.threshold = 1e-5, minlen = 50, minCG = 3,
                dis.merge = 100, pct.sig = 0.5)   # pct.sig=0.5: >=50% of region CpGs individually significant
```

callDMR merges significant CpGs (threshold-then-define), so its region p is NOT selection-corrected; its strengths are dispersion shrinkage and multi-factor support (`DMLtest.multiFactor` with a model formula).

## bsseq BSmooth

**Goal:** Call DMRs on low/uneven-coverage WGBS by smoothing each sample before testing.

**Approach:** Smooth per sample, compute the smoothed t-statistic, then threshold it into regions.

```r
library(bsseq)

bs_smooth <- BSmooth(bs, BPPARAM = MulticoreParam(4), verbose = TRUE)
keep <- rowSums(getCoverage(bs_smooth) >= 2) == ncol(bs_smooth)   # >=2x in every sample
bs_filt <- bs_smooth[keep, ]

tstat <- BSmooth.tstat(bs_filt, group1 = c('C1', 'C2'), group2 = c('T1', 'T2'),
                       estimate.var = 'same', mc.cores = 4)
dmrs <- dmrFinder(tstat, cutoff = c(-4.6, 4.6))   # t-stat cutoff (Hansen 2012 uses quantile-based cutoffs)
```

`dmrFinder` needs `BSmooth.tstat` output, not the smoothed `BSseq` object directly. BSmooth gives a RANKED DMR list, not a calibrated region FDR - dmrseq (same lab lineage) supersedes it for region inference. Over-smoothing (too-wide bandwidth) washes out focal promoter DMRs.

## DMRcate: The Array-vs-WGBS Fork

**Goal:** Call DMRs by Gaussian-kernel smoothing of the per-CpG statistic, with the correct kernel for the platform.

**Approach:** Annotate per-CpG statistics through the array OR sequencing entry point, then smooth and extract.

```r
library(DMRcate)

# ARRAY (450K/EPIC): beta/M matrix; arraytype sets the CpG universe.
design <- model.matrix(~ condition)
ann_array <- cpg.annotate('array', m_values, what = 'M', arraytype = 'EPIC',
                          analysis.type = 'differential', design = design, coef = 2)
dmrs_array <- extractRanges(dmrcate(ann_array, lambda = 1000, C = 2))   # array kernel

# WGBS: different entry point AND a much smaller kernel - the array C=2 over-smooths
# dense sequencing CpGs ~25x. Annotate from a count/edgeR-DSS path, then C=50.
ann_seq <- sequencing.annotate(bs, design = design, coef = 2)
dmrs_seq <- extractRanges(dmrcate(ann_seq, C = 50))   # WGBS kernel (Peters 2021)
```

The default `lambda=1000, C=2` is ARRAY-only; applying it to WGBS produces massively over-smoothed, merged, inflated DMRs. `pcutoff='fdr'` returns no DMRs if the upstream limma/DSS yields no significant CpGs.

## metilene and comb-p (Second Callers)

metilene is a standalone C CLI taking one tab table (chrom, pos, per-sample methylation rate); binary segmentation + a 2D-KS test find data-adaptive boundaries. comb-p is a Python CLI that takes a BED of per-CpG p-values from ANY upstream test, estimates the p-value autocorrelation, does a Stouffer-Liptak-Kechris combination of neighbors, groups regions, and applies a one-step Sidak correction.

```bash
metilene -M 1000 -m 10 -d 0.1 -a g1 -b g2 input.tsv | metilene_output.pl   # -m min CpGs, -d min mean diff
comb-p pipeline -c 4 --seed 0.01 --dist 500 --step 50 -p out methyl_pvals.bed   # --seed = p to start a region
```

Both combine-and-correct rather than re-select; use them as fast corroboration and intersect with dmrseq.

## Thresholds Are Conventions; Region FDR Is Tool-Specific

Every caller exposes the same coupled knobs under different names: min-CpGs (`minNumRegion`/`minCG`/`-m`/`min.cpgs`/`cov.bases`), max-gap (`maxGap`/`dis.merge`/`-M`), delta-beta (`cutoff`/`delta`/`-d`/`betacutoff`/`difference`), and a significance cutoff. Shrinking max-gap, raising min-CpGs, and raising delta all reduce the DMR count, and the same data yields wildly different DMR sets. The phrase "region-level FDR" means three different objects: a selection-aware permutation FDR (dmrseq), a BH/Sidak correction on combined per-CpG p-values (DMRcate/comb-p), or per-unit q on independent tiles/CpGs (methylKit/DSS). Report all knobs, name the tool, and use cross-tool overlap as the evidence statement.

## DMR-to-Gene Mapping and the CpG-Density Enrichment Bias

**Goal:** Interpret DMRs without inflating enrichment from CpG-rich genes.

**Approach:** Annotate DMRs to features (annotatr returns one row per DMR-feature overlap; genomation collapses by precedence), then run enrichment with a method that corrects for CpG/probe count.

```r
library(annotatr)
annots <- build_annotations(genome = 'hg38', annotations = c('hg38_basicgenes', 'hg38_cpg_islands'))
dmr_ann <- annotate_regions(regions = sig, annotations = annots, ignore.strand = TRUE)   # one row per overlap

# Enrichment: methylation has a CpG-density bias (CpG-rich genes harbor DMRs by chance),
# so a plain hypergeometric GO test is biased. Use missMethyl goregion (probe/CpG-bias-aware).
# missMethyl::goregion(sig_ranges, all.cpg=..., collection='GO', array.type='EPIC')
```

A single DMR commonly overlaps or sits between several genes; mapping DMR -> gene (nearest TSS vs overlap vs within-X-kb) is a modeling choice that changes the gene list. Hand the corrected enrichment to pathway-analysis/go-enrichment, flagging that methylation input needs a CpG-bias-aware method (missMethyl `gometh`/`goregion`), not a generic hypergeometric test.

## Single-Sample Domain Structure (NOT Differential)

This is a DIFFERENT problem from the focal between-group callers above. The mammalian methylome partitions at MEGABASE scale into Highly Methylated Domains (HMDs, ~80-90%, ordered) and Partially Methylated Domains (PMDs, ~40-70%, disordered, high-variance), and PMDs coincide with late replication, Lamina-Associated Domains, and the Hi-C B-compartment (Lister 2009 *Nature* 462:315; Berman 2012 *Nat Genet* 44:40). Cancer "global hypomethylation" is a DOMAIN phenomenon - focal CpG-island hypermethylation sitting ON a background of megabase PMD hypomethylation - not a focal one. Domain structure is a SINGLE-SAMPLE, structural question answered by SEGMENTERS, not by any between-group DMR caller (methylKit/DSS/dmrseq/DMRcate/metilene/comb-p have no single-sample segmentation mode).

- **MethylSeekR** (Burger 2013 *Nucleic Acids Res* 41:e155) segments one WGBS methylome into UMRs (CpG-rich unmethylated = promoters/CGIs), LMRs (CpG-poor low-methylated ~30% = distal enhancers), and PMDs. Pipeline: `readMethylome()` -> `plotAlphaDistributionOneChr()` (diagnostic: does the sample have PMDs?) -> `segmentPMDs()` (2-state Gaussian HMM, 101-CpG windows) -> `calculateFDRs()` -> `segmentUMRsLMRs(m=0.5, n=..., pmdGRanges=...)`. PMDs MUST be masked before UMR/LMR calling or PMD disorder spawns spurious LMRs.
- **solo-WCGW** (Zhou 2018 *Nat Genet* 50:591) - an isolated CpG in `[A/T]CG[A/T]` context - loses methylation fastest and most monotonically with cell division and is the most sensitive PMD/mitotic-clock readout, detecting PMD hypomethylation even in near-normal tissue. Quantify as the mean over the published common-PMD solo-WCGW CpG set, not as a DMR. See epigenetic-clocks for the broader clock taxonomy.

The warning: running a focal DMR caller on a PMD-bearing genome manufactures thousands of fake hypo-DMRs that are really one phenomenon - the PMD background shifting - chopped into pieces by the max-gap/min-CpG knobs. Segment domains FIRST, then EXCLUDE PMD intervals from focal calling or STRATIFY every DMR by in-PMD vs out-of-PMD and report the fraction that is PMD background.

## Per-Method Failure Modes

### PMD background reported as DMRs
**Trigger:** focal caller on tumor/aged/placenta/cultured WGBS without domain screening. **Mechanism:** megabase PMD hypomethylation chopped into pieces by max-gap/min-CpG. **Symptom:** thousands of large hypo-DMRs in gene-desert, late-replicating, low-CpG-density coordinates. **Fix:** segment PMDs (MethylSeekR) first; exclude or stratify; report the PMD fraction.

### Pre-smoothing before dmrseq
**Trigger:** `BSmooth()` then feeding the smoothed object to `dmrseq`. **Mechanism:** dmrseq smooths the difference internally; pre-smoothing double-smooths. **Symptom:** distorted candidate regions and invalid statistics. **Fix:** feed dmrseq the raw `BSseq` counts.

### DMRcate array defaults on WGBS
**Trigger:** copying `lambda=1000, C=2` onto sequencing data. **Mechanism:** the array kernel is ~25x too wide for dense WGBS CpGs. **Symptom:** massively over-smoothed, merged, inflated DMRs. **Fix:** `sequencing.annotate()` + `C=50` for WGBS (Peters 2021).

### Threshold-then-test reported as region FDR
**Trigger:** greping runs of significant per-CpG calls and reporting the per-CpG q. **Mechanism:** the regions were selected for extremeness, then tested on the same data. **Symptom:** anti-conservative, uncalibrated region q. **Fix:** use dmrseq (selection-aware) for the headline; at minimum state that a tile/merge q is a screening q.

### Single-CpG tiles
**Trigger:** `tileMethylCounts` at the default `cov.bases=0`. **Mechanism:** a window with one covered CpG becomes a "DMR." **Symptom:** thousands of single-CpG noisy regions. **Fix:** set `cov.bases >= 3`.

### Cross-tool count comparison
**Trigger:** comparing "N DMRs at q<0.05" between callers. **Mechanism:** each tool's q controls a different object. **Symptom:** apparent disagreement that is really an FDR-definition mismatch. **Fix:** compare OVERLAP, not counts.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| dmrseq `cutoff` 0.1 (candidate seed only) | Korthauer 2019 | seeds detection; significance is the permutation statistic, NOT a delta floor |
| DSS `callDMR(delta=0)` default -> SET it | Park & Wu 2016; DSS docs | delta=0 calls regions with no effect-size floor; set ~0.1 to require a real shift |
| min-CpGs per region 3-5 | convention | single-CpG "regions" are DMPs in disguise; trades sensitivity vs specificity |
| delta-beta 25% ("moderate") | methylKit tutorial folklore | NOT derived; biologically meaningful delta is feature- and purity-dependent |
| methylKit `tileMethylCounts(cov.bases>=3)` | nuance (default is 0) | the default 0 lets single-CpG tiles through |
| coverage floor ~10x per CpG | field standard | a single-CpG beta below ~10x is a coin flip |
| DMRcate WGBS `C=50` (array `C=2`) | Peters 2021 *Nucleic Acids Res* 49:e109 | dense WGBS CpGs need a far narrower kernel than array probes |
| MethylSeekR `m=0.5`, FDR<5% | Burger 2013 | methylation cutoff for hypomethylated regions; FDR target picks the CpG-count threshold n |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| dmrseq error about zero-coverage loci | a locus has 0 coverage in some sample | filter `rowSums(getCoverage(bs)==0)==0` first |
| Thousands of huge hypo-DMRs in cancer WGBS | PMD background not segmented | MethylSeekR `segmentPMDs` first; exclude/stratify |
| Over-merged WGBS DMRs with DMRcate | array kernel on sequencing | `sequencing.annotate()` + `C=50` |
| `dmrFinder` errors on a smoothed object | needs `BSmooth.tstat` output | run `BSmooth.tstat` before `dmrFinder` |
| DMRcate returns no DMRs | `pcutoff='fdr'` and no significant upstream CpGs | check the upstream limma/DSS result first |
| Biased GO enrichment of DMR genes | plain hypergeometric ignores CpG density | use missMethyl `goregion`/`gometh` |

## References

- Korthauer K, Chakraborty S, Benjamini Y, Irizarry RA. 2019. Detection and accurate false discovery rate control of differentially methylated regions from whole genome bisulfite sequencing. *Biostatistics* 20:367-383.
- Feng H, Conneely KN, Wu H. 2014. A Bayesian hierarchical model to detect differentially methylated loci from single nucleotide resolution sequencing data. *Nucleic Acids Res* 42:e69.
- Park Y, Wu H. 2016. Differential methylation analysis for BS-seq data under general experimental design. *Bioinformatics* 32:1446-1453.
- Hansen KD, Langmead B, Irizarry RA. 2012. BSmooth: from whole genome bisulfite sequencing reads to differentially methylated regions. *Genome Biol* 13:R83.
- Akalin A, Kormaksson M, Li S, et al. 2012. methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles. *Genome Biol* 13:R87.
- Peters TJ, Buckley MJ, Statham AL, et al. 2015. De novo identification of differentially methylated regions in the human genome. *Epigenetics Chromatin* 8:6.
- Peters TJ, Buckley MJ, Chen Y, et al. 2021. Calling differentially methylated regions from whole genome bisulphite sequencing with DMRcate. *Nucleic Acids Res* 49:e109.
- Juhling F, Kretzmer H, Bernhart SH, Otto C, Stadler PF, Hoffmann S. 2016. metilene: fast and sensitive calling of differentially methylated regions from bisulfite sequencing data. *Genome Res* 26:256-262.
- Pedersen BS, Schwartz DA, Yang IV, Kechris KJ. 2012. Comb-p: software for combining, analyzing, grouping and correcting spatially correlated P-values. *Bioinformatics* 28:2986-2988.
- Lister R, Pelizzola M, Dowen RH, et al. 2009. Human DNA methylomes at base resolution show widespread epigenomic differences. *Nature* 462:315-322.
- Berman BP, Weisenberger DJ, Aman JF, et al. 2012. Regions of focal DNA hypermethylation and long-range hypomethylation in colorectal cancer coincide with nuclear lamina-associated domains. *Nat Genet* 44:40-46.
- Zhou W, Dinh HQ, Ramjan Z, et al. 2018. DNA methylation loss in late-replicating domains is linked to mitotic cell division. *Nat Genet* 50:591-602.
- Burger L, Gaidatzis D, Schubeler D, Stadler MB. 2013. Identification of active regulatory regions from DNA methylation data. *Nucleic Acids Res* 41:e155.

## Related Skills

- differential-cpg-testing - Per-site testing before region aggregation
- methylkit-analysis - methylKit object model and tile construction
- methylation-calling - Produces the input count tables
- array-preprocessing - Array beta/M-value input for DMRcate array mode
- epigenetic-clocks - Mitotic-clock / solo-WCGW overlap (domain section)
- pathway-analysis/go-enrichment - CpG-bias-aware enrichment of DMR genes (missMethyl gometh)
- long-read-sequencing/nanopore-methylation - Pipe modkit bedMethyl counts here for region statistics
- workflows/methylation-pipeline - End-to-end bisulfite pipeline
