# Methylation-Based Detection - Usage Guide

## Overview
Detect cancer and infer tissue-of-origin from cfDNA methylation: choose a conversion chemistry that survives low-input plasma, call read-level methylation haplotypes rather than averaged beta values, and deconvolve a hematopoietic-dominated mixture against a methylation atlas.

## Prerequisites
```bash
conda install -c bioconda methyldackel bismark bwameth
pip install pandas numpy scipy statsmodels
# cfMeDIP enrichment analysis (R/Bioconductor):
# BiocManager::install(c('MEDIPS', 'qsea'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Pick a conversion chemistry for my low-input plasma methylation assay"
- "Run mbias then extract per-CpG methylation from my cfDNA BAM with MethylDackel"
- "Deconvolve tissue-of-origin from cfDNA methylation against a reference atlas"
- "Find region-level DMRs separating cancer from normal cfDNA"
- "Count read-level methylation haplotypes for MRD detection"

## Example Prompts

### Chemistry Choice
> "I have a few nanograms of cfDNA and need base resolution. Should I use bisulfite, EM-seq, or TAPS?"

> "Recommend a genome-wide methylation assay for ultra-low-input plasma and explain the resolution tradeoff."

### Extraction
> "Run MethylDackel mbias to choose trimming, then extract per-CpG methylation with mergeContext from my bisulfite BAM."

> "Parse the MethylDackel CpG bedGraph and compute beta values, accounting for the fixed column order."

### Tissue Deconvolution
> "Deconvolve tissue-of-origin from my cfDNA methylation against the Loyfer atlas using NNLS with the simplex constraint."

> "My tumor coefficient is unstable. Check whether the white-blood-cell background reference is mis-specified."

### DMR Discovery and Detection
> "Find region-level DMRs between cancer and normal cfDNA with Benjamini-Hochberg FDR."

> "Set up read-level methylation haplotype counting over pre-defined blocks for ppm-level MRD detection."

### Enrichment Data
> "Analyze my cfMeDIP-seq data with QSEA to get density-corrected, BS-comparable methylation levels."

## What the Agent Will Do
1. Recommend a conversion chemistry by input, resolution need, and destructiveness
2. Run MethylDackel mbias, then extract per-CpG methylation with the suggested trimming
3. Parse the fixed 6-column bedGraph and compute beta values
4. Discover region-level DMRs with explicit BH FDR (not per-CpG t-tests)
5. Deconvolve tissue-of-origin via NNLS against an atlas with the simplex constraint
6. For detection, count read-level methylation haplotypes rather than averaging beta

## Tips
- **Run mbias first** - feed its suggested --OT/--OB trimming into extract; do not hard-code someone else's bounds
- **Column order is fixed** - bedGraph is chrom / start / end / methylation-% (integer) / count-methylated / count-unmethylated
- **Read-level, not averaged** - one concordant tumor fragment beats a diluted beta mean at low tumor fraction
- **Bisulfite degrades cfDNA** - prefer EM-seq or TAPS for low-input base-resolution work; cfMeDIP for ultra-low input genome-wide
- **cfMeDIP is enrichment** - use MEDIPS/QSEA, never a per-CpG bisulfite pipeline; coverage is not methylation level
- **WBC dominates** - >90% of cfDNA is hematopoietic; the tumor term is a small residual, and a methylation analog of CHIP can confound it
- **Atlas markers are platform-specific** - WGBS-fragment markers do not equal array probes; deconvolve only over markers the assay actually covers

## Related Skills
- cfdna-preprocessing - conversion chemistry and library choices upstream
- fragment-analysis - orthogonal genome-wide cfDNA signal
- analytical-validation - read-level detection framed as a limit-of-detection problem
- methylation-analysis/bismark-alignment - bisulfite read alignment
- methylation-analysis/dmr-detection - region-level differential methylation statistics
