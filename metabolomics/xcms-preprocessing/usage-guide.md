# XCMS Preprocessing - Usage Guide

## Overview

XCMS is the standard Bioconductor package for processing LC-MS metabolomics data. XCMS3 provides a modern, object-oriented interface for peak detection, alignment, and correspondence.

## Prerequisites

```bash
# R/Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("xcms", "MSnbase"))
```

## Quick Start

Tell your AI agent what you want to do:
- "Process my mzML files with XCMS and create a peak table"
- "Detect peaks with CentWave and align samples"

## Example Prompts

### Data Loading
> "Read my mzML files from the data/ folder into an XCMSnExp object"

### Peak Detection
> "Run CentWave peak detection with 10 ppm mass tolerance and 5-20 second peak width"
> "Use MatchedFilter for my profile mode data"

### Alignment
> "Align retention times using Obiwarp with correlation distance function"
> "Apply RT correction across all samples"

### Feature Grouping
> "Group peaks across samples using PeakDensity with 30% minimum sample fraction"
> "Fill in missing peak values after grouping"

### QC Analysis
> "Check the TIC for injection issues and visualize RT alignment"
> "Run PCA on the processed features to check QC sample clustering"

## What the Agent Will Do

1. Load mzML/mzXML files into XCMSnExp object
2. Run peak detection (CentWave or MatchedFilter)
3. Align retention times (Obiwarp)
4. Group corresponding peaks across samples
5. Fill missing values from raw data
6. Export feature table

## Tips

- Use CentWave for centroided data (most modern instruments), MatchedFilter for profile mode
- Set peakwidth based on your chromatography (typically 5-30 seconds for LC)
- Include pooled QC samples every 10 injections for drift correction
- Check TIC plots for injection issues before processing
- QC samples should cluster tightly in PCA

## Key Parameters

| Method | Parameter | Description |
|--------|-----------|-------------|
| CentWave | peakwidth | Expected peak width range (seconds) |
| CentWave | ppm | m/z tolerance in ppm |
| Obiwarp | binSize | m/z bin size for alignment |
| PeakDensity | bw | RT bandwidth for density estimation |
| PeakDensity | minFraction | Min fraction of samples with peak |

## References

- XCMS: doi:10.1021/ac051437y
- XCMS3: doi:10.1021/acs.analchem.7b03003
- Documentation: https://bioconductor.org/packages/xcms/
