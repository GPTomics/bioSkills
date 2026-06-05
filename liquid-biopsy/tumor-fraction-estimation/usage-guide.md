# Tumor Fraction Estimation - Usage Guide

## Overview
Estimate tumor fraction (the proportion of cfDNA molecules that are tumor-derived, the cfDNA analogue of bulk-tumor purity) from shallow whole-genome sequencing with ichorCNA. Tumor fraction is the burden metric that travels across assays; it is not mutation VAF, and CNA-based estimation has a hard ~3 percent floor below which the estimator must change.

## Prerequisites
```r
install.packages('devtools')
devtools::install_github('GavinHaLab/ichorCNA')
```

```bash
# HMMcopy provides readCounter for binning the BAM into WIG
conda install -c bioconda hmmcopy
```

## Quick Start
Tell your AI agent what you want to do:
- "Estimate tumor fraction from my 0.5x sWGS BAM with ichorCNA"
- "Bin my BAM into a 1 Mb WIG and run ichorCNA"
- "My ichorCNA tumor fraction is low but the patient has known disease. Is this a near-diploid false low?"
- "Reconcile my ichorCNA tumor fraction against the max VAF from my panel"
- "Pick a tumor-fraction estimator for a sample I think is below 3 percent"

## Example Prompts

### Single Sample
> "Run readCounter then runIchorCNA.R on my hg38 cfDNA BAM and report tumor fraction and ploidy."

> "Set up the low-tumor-fraction parameter recipe for a sample expected near 1 percent."

### Batch Processing
> "Process all my sWGS WIG files through ichorCNA and collect tumor fractions into one table."

> "Parse the .params.txt across my cohort and flag samples with GC-Map MAD above 0.3."

### Interpretation
> "My ichorCNA tumor fraction is 1.5 percent. Is that real or below the limit of detection?"

> "The genome-wide plot is flat and the tumor fraction is near zero. What does that mean for this near-diploid tumor type?"

> "Convert my panel max VAF to a tumor fraction and check it against ichorCNA."

## What the Agent Will Do
1. Bin the BAM into 1 Mb WIG with readCounter (chromosome naming matched to the BAM @SQ style)
2. Run runIchorCNA.R with build-matched GC/map/centromere references and a protocol-matched panel of normals
3. Estimate tumor fraction (= 1 - n), ploidy, and subclonal prevalence over the normal/ploidy grid
4. Parse .params.txt for the max-loglik solution and the GC-Map MAD QC metric
5. Flag near-diploid false lows, sub-3-percent calls, and ploidy aliasing; recommend an alternative estimator when ichorCNA is out of its regime

## Tips
- Tumor fraction is not VAF: for a clonal heterozygous SNV in a diploid region, tumor fraction is approximately 2 times the VAF.
- The ~3 percent floor is a limit of detection, not a true zero; below it switch to deep-panel max-VAF, methylation deconvolution, or fragmentomics.
- A low tumor fraction is "low burden" only if the genome-wide plot is genuinely flat; in a quiet or copy-neutral-LOH tumor a low value is a false low.
- The panel of normals must match the library prep, bin size, and genome build; a mismatched PoN fabricates copy-number waviness.
- Keep the build (hg19 vs hg38) and chromosome style (1 vs chr1) consistent across the BAM, GC/map/centromere WIGs, and the flags.
- Gate on GC-Map Correction MAD: under 0.15 is good, over 0.3 means distrust the call.
- Read all .params.txt solutions, not just the top one, to catch ploidy-3 aliases of a ploidy-2 truth.

## Related Skills
- cfdna-preprocessing - sWGS BAM input and minimal-processing path
- fragment-analysis - the estimator to use below the ~3% CNA floor
- ctdna-mutation-detection - max-VAF cross-check and the TF-vs-VAF reconciliation
- analytical-validation - the ~3% floor framed as a limit of detection
- copy-number/cnvkit-analysis - copy-number calling concepts
- copy-number/copy-ratio-segmentation - segmentation concepts
