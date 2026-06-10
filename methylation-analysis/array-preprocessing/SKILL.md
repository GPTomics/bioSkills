---
name: bio-methylation-array-preprocessing
description: Turns raw Illumina Infinium methylation BeadChip IDATs (450K, EPIC, EPICv2) into a defensible beta/M matrix with sesame (openSesame/SigDF) or minfi (RGChannelSet -> MethylSet -> GenomicRatioSet). Covers Type I vs Type II probe chemistry and why raw Type II beta is compressed, the signal-to-beta math (beta = M/(M+U+100)) and M-value logit, detection-p / pOOBAH masking including the out-of-band deletion-artifact catch, dye-bias correction, and the normalization decision (noob, funnorm, quantile, SWAN, BMIQ, dasen, sesame QCDPB). Use when reading IDATs, choosing a normalization for a 450K/EPIC/EPICv2 cohort, deciding beta vs M, masking failed probes, or producing the corrected matrix before testing. For probe/sample filtering, EPICv2 replicate collapse, and sample-identity QC see array-qc-filtering; for native long-read 5mC see long-read-sequencing/nanopore-methylation (a different platform).
tool_type: r
primary_tool: sesame
---

## Version Compatibility

Reference examples tested with: sesame 1.20+, minfi 1.48+, ChAMP 2.32+, wateRmelon 2.8+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The ARRAY VERSION and GENOME BUILD are versions that matter as much as the package. EPICv2 REQUIRES sesame (mainstream minfi does not auto-detect it and returns "Unknown"); the manifest/annotation packages are array-version- and genome-build-specific (450K and EPICv1 are hg19; EPICv2 is hg38-native). sesame pulls platform/address data from its own hub, so `sesameDataCache()` must run once before processing. Record the array (450K/EPIC/EPICv2) and the genome build in any output, the way a sequencing run records its reference.

# Array Preprocessing

**"Give me a clean methylation matrix from my IDATs"** -> Read the raw two-channel intensities, correct the Type I/II design mismatch, dye/background bias, and failed probes, then emit beta (for reporting) and M (for testing) - because an Infinium beta is a two-chemistry fluorescence ratio, not a methylation value, until those corrections are applied.
- R: `openSesame(idat_dir, prep='QCDPB', func=getBetas)` (sesame) or `preprocessFunnorm(rgSet)` (minfi)

Scope: IDAT -> corrected, masked, normalized beta/M matrix for one array version. Probe filtering (cross-reactive/SNP/sex), EPICv2 replicate collapse, and sample-identity QC -> array-qc-filtering. Per-CpG testing -> differential-cpg-testing. Region calling -> dmr-detection. Native long-read 5mC -> long-read-sequencing/nanopore-methylation. Bisulfite-sequencing (Bismark/WGBS/RRBS) is the other modality in this category, not this skill.

## The Single Most Important Modern Insight -- An Infinium Beta Value Is a Two-Chemistry Fluorescence Ratio, Not a Methylation Measurement

An Infinium array does not measure methylation - it measures the relative fluorescence of a methylated vs unmethylated allele at a fixed, manufacturer-chosen set of CpGs, glued together from two incompatible chemistries. A raw beta becomes a comparable methylation estimate only after preprocessing; preprocessing IS the measurement, not optional cleanup. Three corollaries every misuse violates:

1. **The manifest is the experiment.** The array interrogates <3% of human CpGs, and a DIFFERENT <3% across 450K (~485K), EPIC (~865K), and EPICv2 (~935K). "Absent" almost always means "not on this array," and a 450K-trained clock or EWAS does not transfer to EPICv2 without intersecting probe sets. Do not start from a supplied beta matrix when IDATs exist - the raw two-channel intensities, control probes, and out-of-band signal that noob/funnorm/pOOBAH need are already gone.
2. **Type I and Type II betas disagree by design.** Type II probes (one bead, two dyes) have a narrower dynamic range and dye-incorporation bias, so raw Type II betas are compressed toward 0.5 relative to Type I (two beads, one channel). Mixing the two chemistries without a design correction (BMIQ/SWAN/sesame matchDesign) injects a probe-type artifact that can exceed the biological effect. The diagnostic is a per-type beta-density plot showing two mismatched peaks.
3. **Raw beta is uninterpretable until detection-masked.** Failed probes - low signal, germline/somatic deletions, cross-reactive, SNP-hit - return confident-looking betas that are pure noise. pOOBAH / detection-p masking is what separates a number from a measurement of nothing; pOOBAH additionally catches deletion-driven false-intermediate methylation that negative-control detection-p misses.

Organize the work around DELIVERING a defensible matrix (read -> correct -> mask -> normalize), not around listing minfi functions.

## Three Modalities of the Same Biology

DNA methylation is measured three ways, each with different tradeoffs - state which one the data is before choosing tools:

| Modality | Readout | Coverage | Cohort-comparability | This skill |
|----------|---------|----------|----------------------|------------|
| Infinium array (450K/EPIC/EPICv2) | intensity ratio, no depth | fixed <3% of CpGs, regulatory-enriched | high (shared manifest, no alignment) | YES |
| WGBS / RRBS bisulfite | count ratio, depth-gated | genome-wide (WGBS) or enriched (RRBS) | needs alignment + matched genome | -> bismark-alignment, methylkit-analysis |
| Long-read native (ONT/PacBio) | per-molecule modification calls | genome-wide, phased | growing | -> long-read-sequencing/nanopore-methylation |

Arrays dominate human epigenetic epidemiology (essentially every published clock and large EWAS is array-based) because cost is a fraction of WGBS and the fixed manifest makes cohorts directly comparable.

## Object Models (do not start from a beta matrix)

The raw output per sample is a pair of binary IDATs (`_Grn.idat`, `_Red.idat`); background, dye, and detection-p correction REQUIRE these plus the control probes.

- **minfi:** `RGChannelSet` (raw red/green) -> a `preprocess*` step -> `MethylSet` (M/U intensities) -> `RatioSet` (beta/M) -> `GenomicRatioSet` (genome-mapped). `read.metharray.exp()` reads IDATs; `getBeta()`, `getM()`, `getCN()` extract values.
- **sesame:** a `SigDF` (one signal data.frame per sample). `readIDATpair()` reads one sample; `openSesame()` drives the whole pipeline across a directory and returns a betas matrix directly.

## Tool Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| sesame | Zhou 2018 *Nucleic Acids Res* 46:e123 | SigDF; openSesame QCDPB; pOOBAH OOB masking; EPICv2-native | EPICv2; best detection masking; the modern default |
| minfi | Aryee 2014 *Bioinformatics* 30:1363 | RGChannelSet->GenomicRatioSet; noob/funnorm/quantile/SWAN | 450K/EPICv1; large downstream ecosystem (DMRcate, conumee) |
| ChAMP | Tian 2017 *Bioinformatics* 33:3982 | end-to-end pipeline; BMIQ default | one-call newcomer pipeline on 450K/EPICv1 |
| wateRmelon | Pidsley 2013 *BMC Genomics* 14:293 | dasen/nasen + metric-driven normalization eval | dasen default; normalization benchmarking |

## Normalization Decision Tree by Scenario

Separate the two correction layers that get conflated: (a) background + dye bias (within-sample): noob, sesame dyeBias, dasen background step; (b) Type I/II design correction + between-array harmonization: SWAN, BMIQ, quantile, funnorm, dasen quantile step. A complete pipeline does both.

| Scenario | Recommended | Why |
|----------|-------------|-----|
| EPICv2 (any design) | sesame `openSesame(prep='QCDPB')` | EPICv2-native; pOOBAH; minfi mis-handles duplicate IDs |
| Cancer / cross-tissue (global differences expected) | minfi `preprocessFunnorm` (noob + control-PCs) | preserves real global shifts; quantile would erase them |
| Subtle blood EWAS (no global difference expected) | `preprocessQuantile` or wateRmelon `dasen` | marginal distributions assumed equal; safe to harmonize |
| Strong Type I/II design correction wanted | BMIQ (Teschendorff 2013) or SWAN (Maksimovic 2012) | dilate Type II onto the Type I distribution; pair with a between-array step |
| Single-sample / clinical / streaming | ssNoob or per-IDAT openSesame | reproducible without re-normalizing the cohort |
| Probe/sample filtering, EPICv2 collapse, identity | -> array-qc-filtering | this skill stops at the corrected matrix |
| Per-CpG testing on the matrix | -> differential-cpg-testing | test on M-values; report delta-beta |

There is no universally best normalization (Pidsley 2013 favored dasen; Fortin 2014 favored funnorm for global-difference studies; Welsh 2023 ranked a sesame/pOOBAH pipeline best and quantile worst on EPIC replicate-concordance). Key the choice on array version + whether global differences are expected + single-sample vs cohort, and verify against current benchmarks rather than hard-coding one method.

## Signal -> Beta -> M

- **Beta:** `beta = M / (M + U + alpha)`, M = methylated-allele intensity, U = unmethylated, `alpha = 100` (minfi default) stabilizes the ratio when both intensities are near zero. beta in [0,1] is interpretable but HETEROSCEDASTIC (variance collapses near 0 and 1), violating the constant-variance assumption of linear models.
- **M-value:** `M = log2((M_int + alpha) / (U_int + alpha))`, the logit of beta. Approximately homoscedastic; the correct scale for limma/t-tests (Du 2010 *BMC Bioinformatics* 11:587). Rule: test on M-values, report delta-beta for effect size - the same rule as bisulfite sequencing.

## Process IDATs with sesame (the EPICv2-safe default)

**Goal:** Produce a corrected, detection-masked betas matrix from a directory of IDAT pairs without manually juggling manifest packages.

**Approach:** Cache the sesame data hub once, then run openSesame with the default `QCDPB` prep (qualityMask, inferInfiniumIChannel, dyeBiasNL, pOOBAH, noob, in that order), which auto-detects the platform and returns betas; pOOBAH writes NA into failed probes in place.

```r
library(sesame)
sesameDataCache()                          # once per machine; pulls platform/address data
betas <- openSesame('idat_dir', prep = 'QCDPB', func = getBetas)
# prep codes: Q qualityMask  C inferInfiniumIChannel  D dyeBiasNL  P pOOBAH  B noob
# pOOBAH masks (sets NA) probes whose out-of-band signal is indistinguishable from background,
# catching deletion-driven false-intermediate methylation that negative-control detection-p misses
mvals <- log2(betas / (1 - betas))         # M-values for statistical testing (logit of beta)
```

For EPICv2, openSesame detects the platform automatically; the replicate-probe collapse (`betasCollapseToPfx`) belongs to the next stage and is documented in array-qc-filtering.

## Process IDATs with minfi (450K / EPICv1)

**Goal:** Build a normalized GenomicRatioSet and extract beta and M, choosing the normalization by whether global methylation differences are expected.

**Approach:** Read IDATs into an RGChannelSet, compute a detection-p mask before normalizing, then apply funnorm (global differences) or quantile (no global differences); extract beta and M with the offset-100 defaults.

```r
library(minfi)
rgSet <- read.metharray.exp(base = 'idat_dir')
detP <- detectionP(rgSet)                   # neg-control-based; pre-normalization probe-failure map

grSet <- preprocessFunnorm(rgSet, nPCs = 2) # noob first, then 2 control-probe PCs; preserves global shifts
# preprocessQuantile(rgSet) instead when NO global difference is expected (subtle blood EWAS)

beta <- getBeta(grSet)                       # GenomicRatioSet holds precomputed betas (offset applied upstream)
mval <- getM(grSet)                          # log2(beta/(1-beta)) on the ratio set
beta[detP[rownames(beta), colnames(beta)] > 0.01] <- NA   # mask probes failing detection-p (0.01)
```

EPICv2 is NOT handled by mainstream minfi (it returns "Unknown" and duplicates probe IDs); use sesame for EPICv2.

## Per-Method Failure Modes

### Starting from a supplied beta matrix
**Trigger:** processing begins from a `.csv`/`.RData` beta matrix instead of IDATs. **Mechanism:** a beta matrix has discarded the raw two-channel intensities, control probes, and out-of-band signal. **Symptom:** noob/funnorm/pOOBAH/dye correction cannot run; detection-p cannot be recomputed. **Fix:** obtain the raw IDAT pairs; treat a beta matrix as a last resort and document that preprocessing could not be applied.

### Type I/II mismatch left in the data
**Trigger:** testing on raw or only background-corrected betas. **Mechanism:** Type II betas are compressed toward 0.5 relative to Type I. **Symptom:** "differential" probes that are design artifacts; a two-peak per-type beta density. **Fix:** apply BMIQ/SWAN or sesame matchDesign (or use openSesame, which corrects channel/dye) before testing; confirm the two per-type peaks align.

### minfi on EPICv2
**Trigger:** `read.metharray.exp` on EPICv2 IDATs with mainstream minfi. **Mechanism:** EPICv2 is not auto-detected; 5,483 loci carry duplicate IDs. **Symptom:** array reads as "Unknown"; `getBeta()` returns repeated rownames so `match()`-based merges silently misbehave. **Fix:** use sesame (EPICv2-native), or install a third-party EPICv2 manifest/anno, tag the annotation manually, and collapse replicates in array-qc-filtering.

### Quantile-normalizing a global-difference contrast
**Trigger:** `preprocessQuantile` on cancer vs normal or cross-tissue data. **Mechanism:** between-array quantile assumes equal marginal beta distributions. **Symptom:** real global hypomethylation flattened away. **Fix:** use funnorm (control-probe PCs preserve global shifts); reserve quantile/dasen for subtle no-global-difference designs.

### Testing on beta instead of M
**Trigger:** limma/t-tests run directly on beta. **Mechanism:** beta is heteroscedastic (variance collapses near 0 and 1). **Symptom:** miscalibrated variance; inflated or deflated p-values at extreme methylation. **Fix:** test on M-values, report delta-beta for effect size.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| beta offset `alpha = 100` | Aryee 2014; minfi default | stabilizes the ratio when M and U are both near zero |
| detection-p `> 0.01` = failed | minfi convention | signal indistinguishable from background; beta is noise |
| pOOBAH default p ~ 0.05 | Zhou 2018 | OOB-based mask; also catches deletion-driven false intermediate methylation |
| funnorm `nPCs = 2` | Fortin 2014 | first 2 control-probe PCs absorb technical variation without erasing biology |
| Test on M-values, report delta-beta | Du 2010 | M is homoscedastic for modeling; beta is interpretable for effect size |
| sesame prep `QCDPB` (ordered) | Zhou 2018 | Q quality, C channel, D dye, P pOOBAH, B noob - the validated default order |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Array reads as "Unknown" | EPICv2 in mainstream minfi | use sesame; or third-party manifest + manual annotation tag |
| `getBeta()` has repeated rownames | EPICv2 duplicate probe IDs | collapse replicates (array-qc-filtering); do not merge by ID first |
| Two-peak beta density per probe type | Type I/II design bias uncorrected | BMIQ/SWAN/openSesame before testing |
| Global signal vanished after normalization | quantile applied to a global-difference study | use funnorm |
| `sesameDataCache` / platform-not-found | hub not cached | run `sesameDataCache()` once before processing |
| Coordinates misalign merging EPICv2 with 450K | EPICv2 is hg38, 450K/EPICv1 hg19 | track build per array; liftover before merging (array-qc-filtering) |

## References

- Aryee MJ, Jaffe AE, Corrada-Bravo H, et al. 2014. Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays. *Bioinformatics* 30:1363-1369.
- Zhou W, Triche TJ Jr, Laird PW, Shen H. 2018. SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions. *Nucleic Acids Res* 46:e123.
- Triche TJ Jr, Weisenberger DJ, Van Den Berg D, Laird PW, Siegmund KD. 2013. Low-level processing of Illumina Infinium DNA methylation BeadArrays. *Nucleic Acids Res* 41:e90.
- Fortin JP, Labbe A, Lemire M, et al. 2014. Functional normalization of 450k methylation array data improves replication in large cancer studies. *Genome Biol* 15:503.
- Maksimovic J, Gordon L, Oshlack A. 2012. SWAN: subset-quantile within array normalization for Illumina Infinium HumanMethylation450 BeadChips. *Genome Biol* 13:R44.
- Teschendorff AE, Marabita F, Lechner M, et al. 2013. A beta-mixture quantile normalization method for correcting probe design bias in Illumina Infinium 450k DNA methylation data. *Bioinformatics* 29:189-196.
- Pidsley R, Wong CCY, Volta M, Lunnon K, Mill J, Schalkwyk LC. 2013. A data-driven approach to preprocessing Illumina 450K methylation array data. *BMC Genomics* 14:293.
- Du P, Zhang X, Huang CC, et al. 2010. Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis. *BMC Bioinformatics* 11:587.
- Tian Y, Morris TJ, Webster AP, et al. 2017. ChAMP: updated methylation analysis pipeline for Illumina BeadChips. *Bioinformatics* 33:3982-3984.
- Kaur D, Lee SM, Goldberg D, et al. 2023. Comprehensive evaluation of the Infinium human MethylationEPIC v2 BeadChip. *Epigenetics Commun* 3:6.

## Related Skills

- array-qc-filtering - Probe and sample QC/filtering downstream of preprocessing
- differential-cpg-testing - Per-CpG testing on the resulting beta/M matrix
- dmr-detection - DMRcate array-mode region calling
- cell-type-deconvolution - Consumes the clean beta matrix
- epigenetic-clocks - Consumes the clean beta matrix
- ewas-design - Study design, batch, and inference layer
- long-read-sequencing/nanopore-methylation - Native long-read methylation (different platform)
- workflows/methylation-pipeline - End-to-end pipeline
