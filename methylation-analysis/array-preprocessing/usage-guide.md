# Array Preprocessing

## Overview

Illumina Infinium methylation BeadChips (450K, EPIC, EPICv2) are the dominant platform for human epigenetic epidemiology, but a raw beta value is a two-chemistry fluorescence ratio conditioned on a fixed manifest, not a methylation measurement. This skill guides an AI agent through turning raw IDAT files into a defensible beta/M matrix: reading the two-channel intensities, correcting the Type I vs Type II probe-design mismatch, removing dye and background bias, masking failed probes (detection-p / pOOBAH), and choosing a normalization keyed on the array version and whether global methylation differences are expected. It stops at the corrected matrix; probe filtering, EPICv2 replicate collapse, and sample-identity QC belong to array-qc-filtering.

## Prerequisites

R with Bioconductor and the array-version-specific packages:

```r
BiocManager::install(c('sesame', 'sesameData', 'minfi', 'ChAMP', 'wateRmelon'))
# 450K:    IlluminaHumanMethylation450kmanifest, IlluminaHumanMethylation450kanno.ilmn12.hg19
# EPICv1:  IlluminaHumanMethylationEPICmanifest, IlluminaHumanMethylationEPICanno.ilm10b4.hg19
# EPICv2:  requires sesame (mainstream minfi does not auto-detect it)
BiocManager::install(c('IlluminaHumanMethylation450kmanifest',
                       'IlluminaHumanMethylation450kanno.ilmn12.hg19',
                       'IlluminaHumanMethylationEPICmanifest',
                       'IlluminaHumanMethylationEPICanno.ilm10b4.hg19'))
```

Conceptual prerequisites:

- Raw IDAT pairs (`*_Grn.idat`, `*_Red.idat`) plus a sample sheet. Do not start from a supplied beta matrix when IDATs exist; correction needs the raw intensities and control probes.
- `sesameDataCache()` must run once before sesame processing to pull platform/address data.
- The array version (450K / EPIC / EPICv2) and genome build (hg19 for 450K/EPICv1, hg38 for EPICv2) must be known and recorded.
- EPICv2 requires sesame; minfi is fine for 450K and EPICv1.

## Quick Start

Tell your AI agent what you want to do:

- "Read my EPICv2 IDATs and give me a corrected beta and M-value matrix"
- "Preprocess this 450K cancer cohort and preserve the global hypomethylation"
- "Mask failed probes with pOOBAH and return the betas"
- "Which normalization should I use for a subtle blood EWAS?"

## Example Prompts

### EPICv2 preprocessing

> "I have a directory of EPICv2 IDAT pairs. Cache the sesame data, run openSesame with the default QCDPB prep to mask failed probes and correct dye and channel bias, and return both a beta matrix and an M-value matrix."

### Choosing a normalization

> "My study compares tumor and adjacent-normal tissue on 450K arrays. I expect global hypomethylation. Read the IDATs with minfi, normalize in a way that does not erase the global difference, and extract beta and M."

### Detection masking

> "Process these EPIC IDATs, compute a detection p-value per probe, and set probes failing at p > 0.01 to NA before you hand me the matrix."

### Diagnosing the Type I/II artifact

> "After preprocessing, plot the beta density split by Infinium design type to confirm the Type I and Type II peaks now line up, and tell me whether a BMIQ step is still needed."

## What the Agent Will Do

1. Confirm the array version and genome build, and that raw IDAT pairs (not a beta matrix) are available.
2. Choose the framework: sesame for EPICv2 or best detection masking; minfi for 450K/EPICv1 when the downstream ecosystem is needed.
3. Read the IDATs into the raw object (sesame SigDF via openSesame, or minfi RGChannelSet via `read.metharray.exp`).
4. Apply correction: channel inference, dye-bias, background (noob), and detection masking (pOOBAH or detection-p at 0.01).
5. Select a normalization keyed on the scenario: funnorm for expected global differences, quantile/dasen for subtle no-global-difference designs, BMIQ/SWAN for strong Type I/II correction, ssNoob/openSesame for single-sample.
6. Extract beta (for reporting, with offset 100) and M-values (logit, for testing).
7. Record the array version, genome build, normalization method, and probe-masking outcome; hand off filtering and identity QC to array-qc-filtering.

## Tips

- Always process from IDATs. A supplied beta matrix has thrown away the control probes and out-of-band signal that every correction depends on.
- Test on M-values, report delta-beta. Beta is heteroscedastic; M is the homoscedastic scale for linear models.
- The funnorm-vs-quantile choice hinges on one question: are global methylation differences expected? Cancer and cross-tissue contrasts say yes (use funnorm); a subtle blood exposure says no (quantile/dasen is safe).
- For EPICv2, use sesame. Mainstream minfi returns "Unknown" and duplicates probe IDs, which silently breaks ID-based merges.
- A two-peak per-type beta density after preprocessing means the Type I/II design bias is still present; add BMIQ or SWAN (or use openSesame, which corrects channel and dye).
- EPICv2 is hg38; 450K and EPICv1 are hg19. Track the build and liftover before merging coordinates across versions.

## Related Skills

- array-qc-filtering - Probe and sample QC/filtering downstream of preprocessing
- differential-cpg-testing - Per-CpG testing on the resulting beta/M matrix
- dmr-detection - DMRcate array-mode region calling
- cell-type-deconvolution - Consumes the clean beta matrix
- epigenetic-clocks - Consumes the clean beta matrix
- ewas-design - Study design, batch, and inference layer
- long-read-sequencing/nanopore-methylation - Native long-read methylation (different platform)
- workflows/methylation-pipeline - End-to-end pipeline
