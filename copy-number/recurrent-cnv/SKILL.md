---
name: bio-copy-number-recurrent-cnv
description: Identify recurrent and driver copy number alterations across a tumor cohort with GISTIC2 (G-score, Ziggurat deconstruction, focal vs broad/arm-level analysis, q-values from permutation) and quantify copy-number signatures with the Steele 2022 COSMIC framework and the Drews 2022 CINSignatures framework. Covers driver-gene localization from recurrence peaks, distinguishing focal drivers from arm-level passengers, and the caller-sensitivity caveats of copy-number signatures. Use when finding recurrently amplified or deleted regions in a cohort, localizing driver genes, separating focal from broad events, running GISTIC2, or extracting copy-number mutational signatures.
tool_type: mixed
primary_tool: gistic2
---

## Version Compatibility

Reference examples tested with: GISTIC 2.0.23, R 4.3+ with CINSignatureQuantification 1.2+; Python 3.10+ with SigProfilerAssignment 0.1+ (optional, COSMIC CN signatures).

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `gistic2 --help` (GISTIC 2.0 is a MATLAB-compiled binary; needs the MCR runtime)
- R: `packageVersion('CINSignatureQuantification')`
- Python: `pip show SigProfilerAssignment`

GISTIC 2.0 has had no substantive release since ~2017; it is effectively frozen. It runs as a compiled binary against the MATLAB Compiler Runtime — there is no R or Python package. Verify the reference (`-refgene`) `.mat` file matches the genome build.

# Recurrent and Driver Copy Number Alteration

**"Which copy number changes recur across my cohort, and which gene is the driver"** -> A CNV in one tumor is an observation; a CNV recurring across many tumors beyond chance is evidence of selection. GISTIC2 separates recurrent driver events from passengers by modeling a background rate and scoring each locus by how often, and how strongly, it is altered. Copy-number signatures decompose the genome-wide pattern of alterations into the mutational processes that generated them.

- CLI: `gistic2` — cohort-level recurrence, focal vs broad, driver localization
- R: `CINSignatureQuantification` (Drews 2022); Python `SigProfilerAssignment` (Steele 2022 COSMIC)

## How GISTIC2 Works — and Its Limits

GISTIC2 scores each genomic marker with a **G-score** = frequency of alteration x mean amplitude, separately for amplifications and deletions. Significance (**q-value**) comes from permuting events along the genome under the null that all are passengers. **Ziggurat deconstruction** decomposes each sample's profile into the additive arm-level and focal events that produced it, so the background rate is estimated separately for broad and focal alterations — without this, ubiquitous arm-level events swamp the focal signal. A **peel-off** procedure removes the contribution of each significant peak before testing the next, so one strong driver does not mask its neighbors.

Two postdoc-level caveats define how GISTIC2 output must be read:

1. **q-values are cohort-size dependent.** Larger N manufactures more "significant" peaks. A peak list from N=50 and one from N=500 are not comparable; recurrence *frequency* is the portable quantity, not the q-value.
2. **GISTIC2 is only as good as its input segmentation.** Oversegmented seg files produce spurious narrow peaks. The seg file must also be correctly **centered** on diploid — a mis-centered profile (WGD genome centered on tetraploid) inverts every call before GISTIC even runs.

## Decision Tree

| Goal | Approach | Notes |
|------|----------|-------|
| Find recurrent focal drivers in a cohort | GISTIC2, focal analysis, peak regions | Driver = recurrence-peak gene with a known role |
| Quantify arm-level / broad events | GISTIC2 `-broad 1`, arm-level output | `-brlen` sets the focal/broad length cutoff |
| Compare cohorts of different size | Recurrence frequency, not q-value | q-value is not portable across N |
| Characterize mutational processes | Copy-number signatures | Drews CINSignatures or Steele COSMIC CN |
| Localize the gene within a wide peak | GISTIC2 `-genegistic 1` + known drivers | Wide peaks need orthogonal driver evidence |
| Single tumor (no cohort) | GISTIC2 does not apply | Use focal-amplification-ecdna / per-sample annotation |

## Running GISTIC2

```bash
# Segment file: 6 columns -- sample, chrom, start, end, num_markers, seg.mean (log2).
# It MUST be diploid-centered. Pool per-sample segments (e.g. cnvkit.py export seg).
gistic2 \
    -b gistic_output/ \
    -seg cohort.seg \
    -refgene hg38.refgene.mat \
    -genegistic 1 \
    -broad 1 \
    -brlen 0.7 \
    -conf 0.99 \
    -armpeel 1 \
    -savegene 1 \
    -gcm extreme \
    -rx 0
```

Key flags: `-brlen 0.7` sets the focal/broad cutoff at 70% of a chromosome arm; `-conf 0.99` is the peak-boundary confidence — raising it above the 0.75 default yields a wider, more conservative peak with higher confidence the true driver gene lies inside it (the trade-off is more genes per peak); `-armpeel 1` peels arm-level events before focal testing; `-genegistic 1` runs the gene-level test; `-rx 0` keeps sex chromosomes. Output `amp_genes.txt` / `del_genes.txt` and `all_lesions.txt` list peaks, q-values, and genes.

## Copy-Number Signatures

**Goal:** Decompose the genome-wide copy-number pattern into mutational processes (HRD, chromothripsis, tandem duplication, ecDNA, whole-genome doubling).

**Approach:** Two competing 2022 frameworks exist. Steele et al (Nature 2022) defined 21 pan-cancer CN signatures from a 48-channel feature matrix, now in COSMIC; Drews et al (Nature 2022) defined 17 signatures via the CINSignatures feature set. Quantify against one framework consistently; signatures require *absolute* (allele-specific) copy number.

```r
library(CINSignatureQuantification)

# segments: data frame with columns chromosome, start, end, segVal (total CN),
# sample -- absolute copy number from ASCAT/Sequenza/FACETS, NOT relative log2.
res <- quantifyCNSignatures(segments, experimentName = 'cohort',
                            method = 'drews')
activities <- getActivities(res)   # samples x signatures exposure matrix
```

The critical caveat (Steele 2022): three signatures had to be discarded as oversegmentation artifacts and ten were linear combinations needing manual filtering. Signatures are sensitive to the upstream caller — Steele prescribes the caller per platform (SNP6 -> ASCAT penalty 70; shallow WGS -> ASCAT.sc) precisely for this reason.

## Failure Modes

### Comparing q-values across cohorts of different size

**Trigger:** Stating that cohort A has "more significant" peaks than cohort B when the cohorts differ in N.

**Mechanism:** GISTIC q-values fall as N rises — the same recurrence frequency clears significance in a larger cohort.

**Symptom:** A larger cohort appears to have more drivers purely because it is larger; peak lists do not replicate.

**Fix:** Compare recurrence *frequency* (fraction of samples altered), not q-value, across cohorts. Re-run GISTIC at matched N (subsample) if a significance comparison is unavoidable.

### Oversegmented input produces spurious peaks

**Trigger:** Feeding GISTIC a seg file from a noisy or over-fragmented segmentation.

**Mechanism:** GISTIC interprets every segment edge as a potential focal event boundary; fragmentation creates many narrow false peaks.

**Symptom:** Numerous tiny significant peaks at no known driver; peaks not replicated with a cleaner segmentation.

**Fix:** Quality-control the segmentation first (see copy-ratio-segmentation); merge over-fragmented segments before pooling the cohort seg file.

### Mis-centered seg file inverts everything

**Trigger:** Pooling seg files that are not diploid-centered (e.g. WGD tumors centered on tetraploid).

**Mechanism:** GISTIC assumes seg.mean ~ 0 is diploid; a shifted baseline turns gains into neutral and neutral into losses before any statistics run.

**Symptom:** Amplification and deletion peaks swapped relative to known biology; genome-wide deletion bias.

**Fix:** Center each sample's seg file on its true diploid baseline (anchor with allele-specific ploidy) before pooling. Do not rely on per-sample median centering for aneuploid cohorts.

### Treating a wide GISTIC peak as a single-gene call

**Trigger:** Reporting every gene inside a wide significant peak, or assuming the peak gene is the driver.

**Mechanism:** Peak width reflects breakpoint heterogeneity across the cohort; a wide peak may contain dozens of genes, and the statistical peak need not coincide with the functional driver.

**Symptom:** A multi-gene peak reported as one driver; the named gene is a passenger.

**Fix:** Intersect peaks with known drivers (COSMIC CGC, OncoKB), expression, and dependency data. Raising `-conf` widens the peak (it does not narrow it) — peak width is set by cohort breakpoint heterogeneity, not a tunable. Wide peaks require orthogonal driver evidence — GISTIC localizes, it does not nominate.

### Copy-number signatures from relative copy number

**Trigger:** Running CN signatures on log2 ratios or relative segments.

**Mechanism:** Signature features (segment size, copy-number state, change-point) are defined on absolute copy number; relative input gives meaningless states.

**Symptom:** Implausible signature exposures; ploidy/WGD signatures fire spuriously.

**Fix:** Use absolute allele-specific copy number from ASCAT/Sequenza/FACETS as input. Apply the framework's prescribed caller for the platform.

## Reconciliation

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| GISTIC peak with no known driver | Wide peak, passenger locus, or fragile site | Cross-check expression/dependency; treat as candidate |
| Focal peak inside a broad event | Arm-level event not peeled | Confirm `-armpeel 1`; inspect Ziggurat output |
| Drews vs Steele signatures disagree | Different feature definitions and reference sets | Pick one framework; do not mix exposures |
| Peaks change with segmentation | Input over/under-segmented | Stabilize segmentation; re-run |

**Operational rule:** Report a GISTIC peak as a candidate driver locus only when (1) the input segmentation is QC-passed and diploid-centered, (2) recurrence frequency (not just q-value) is substantial, and (3) the peak contains a gene with independent driver evidence. Signatures are reportable only from absolute CN with a single, platform-matched framework.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| GISTIC significance | q < 0.25 | GISTIC2 default residual-q cutoff for peaks |
| Peak-boundary confidence | `-conf 0.99` | Wider, conservative peak; higher confidence the true driver is inside (default 0.75) |
| Focal/broad cutoff | `-brlen 0.7` | Events > 70% of an arm are treated as broad |
| Cohort size for stable peaks | tens to hundreds | Mehta-style: too few samples gives unstable peaks; q is N-dependent |
| CN signatures input | absolute (allele-specific) CN | Steele 2022 / Drews 2022; relative log2 is invalid |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| GISTIC2 will not start | MATLAB Compiler Runtime missing | Install the MCR version GISTIC was built against |
| Amp/del peaks swapped vs biology | Seg file not diploid-centered | Center on true ploidy before pooling |
| Many tiny spurious peaks | Oversegmented input | QC and merge segmentation first |
| `-refgene` errors | Build mismatch (hg19 vs hg38 .mat) | Use the matching reference .mat |
| Implausible signature exposures | Relative CN used as input | Use absolute allele-specific CN |
| Peak lists do not replicate | q-value compared across different N | Compare recurrence frequency |

## References

- Mermel CH et al 2011. GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers. Genome Biol 12:R41
- Beroukhim R et al 2010. The landscape of somatic copy-number alteration across human cancers. Nature 463:899
- Steele CD et al 2022. Signatures of copy number alterations in human cancer. Nature 606:984
- Drews RM et al 2022. A pan-cancer compendium of chromosomal instability. Nature 606:976
- Macintyre G et al 2018. Copy number signatures and mutational processes in ovarian carcinoma. Nat Genet 50:1262

## Related Skills

- copy-number/allele-specific-copy-number - Absolute CN input for GISTIC and CN signatures
- copy-number/copy-ratio-segmentation - Segmentation quality controlling GISTIC peaks
- copy-number/cnv-annotation - Annotating GISTIC peaks with genes and driver roles
- copy-number/focal-amplification-ecdna - Resolving the architecture of focal amplicons
- copy-number/cnv-visualization - Cohort heatmaps of recurrent CNV
- pathway-analysis/go-enrichment - Pathway context for recurrently altered genes
