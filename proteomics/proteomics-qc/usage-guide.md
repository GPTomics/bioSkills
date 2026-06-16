# Proteomics QC - Usage Guide

## Overview
Quality control for bottom-up proteomics framed as a three-level funnel: instrument/raw-signal QC (mass accuracy, RT/iRT fit, FWHM, TIC vs injection time, % MS2 identified), identification/run QC (missed cleavages, charge states, sample-handling PTM artifacts, contaminants), and experiment/quantitative QC (replicate correlation on log2, CV on the linear scale, completeness, missingness mechanism, PCA/batch, TMT channel balance, DIA q-values). The deliverable matrix is the LAST place a fault becomes visible, so the agent inspects raw signal and removes contaminants BEFORE normalizing, because median normalization erases loading failures.

## Prerequisites
```bash
pip install numpy pandas scipy scikit-learn matplotlib seaborn
# R packages: install.packages('PTXQC')   # PTXQC is on CRAN, NOT Bioconductor
# BiocManager::install(c('limma', 'MSstatsTMT'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Plot raw per-sample boxplots and ID counts before I normalize, and flag loading failures"
- "Strip MaxQuant contaminant and decoy rows before log-transform and normalization"
- "Compute replicate correlation on log2 and CV on the linear scale, then run PCA colored by batch"
- "Diagnose whether my missing values are MNAR or MCAR before I pick an imputer"

## Example Prompts

### Inspect Before Normalizing
> "Show raw, un-normalized boxplots per sample with ID counts and total signal, and flag any sample shifted more than 2-3x below its group"

> "Remove rows flagged Potential contaminant, Reverse, or Only identified by site from my MaxQuant proteinGroups before I normalize"

> "What fraction of summed intensity is keratin and trypsin, and is it higher in my low-input samples?"

### Reproducibility
> "Calculate within-group Pearson correlation on log2 intensities and tell me if any sample correlates better with a different group"

> "Compute median CV per condition on the linear scale, not on log data"

> "Are my technical replicates above r 0.98 and is the biological CV in the 20-40% range?"

### Missing Values
> "Plot present-fraction versus abundance to decide if missingness is MNAR or MCAR"

> "Filter to proteins valid in at least 70% of replicates in one condition before imputing"

> "Which imputation method matches my missingness mechanism, and why would the wrong one corrupt my results?"

### Batch and Outliers
> "Run PCA on the normalized survivors and test whether PC1 or PC2 associates with batch rather than condition"

> "Is batch the dominant axis of variance, and should I correct it before differential testing?"

### TMT and DIA
> "Check TMT channel-loading balance on raw reporter intensities and flag any channel deviating more than 2x"

> "How many protein groups pass at 1% global q-value, and why is precursor q not enough?"

## What the Agent Will Do
1. Load the un-normalized search output and intensity matrix
2. Plot raw per-sample boxplots, ID counts, total signal, and missing fraction; flag loading/injection failures
3. Filter contaminant, reverse, and only-identified-by-site rows before any transform
4. Log2-transform and normalize on the surviving good samples (mechanics route to quantification)
5. Compute replicate Pearson r on log2 and median CV on the linear scale
6. Diagnose the missingness mechanism and apply a completeness filter before imputing
7. Run PCA, test for batch association, and flag outlier samples
8. For TMT, inspect channel balance and isolation interference; for DIA, summarize global q-value counts

## Tips
- Inspect raw signal first; once the matrix is median-normalized the evidence of a loading failure is mathematically gone.
- Almost no metric has a universal cutoff; trend each metric against your own per-instrument rolling baseline (Levey-Jennings, +/-2 SD warn, +/-3 SD action).
- Read every metric with its co-readouts: a protein-count drop means spray, column, or sample depending on what co-moves.
- Correlate on log2, never raw; a few abundant proteins make raw r meaningless.
- Compute CV on linear intensity (or geometric CV on logs); always state normalization, transform, and software or the CV is uninterpretable.
- Diagnose MNAR vs MCAR before choosing an imputer; the wrong one either kills or fabricates differences.
- Document and justify every sample exclusion, and run a sensitivity check with and without borderline samples.

## Related Skills
- data-import - Load search-engine output and intensity matrices before QC
- quantification - Normalization and imputation mechanics that QC mandates running AFTER inspection
- differential-abundance - The moderated statistical test QC gates
- dia-analysis - DIA q-value/FDR internals behind the protein-count QC
- data-visualization/dimensionality-reduction-plots - PCA/MDS projection plotting
- workflows/proteomics-pipeline - End-to-end pipeline placing QC before differential testing
