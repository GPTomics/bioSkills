# XCMS Preprocessing Usage Guide

## Overview

XCMS is the standard Bioconductor package for processing LC-MS metabolomics data. XCMS3 provides a modern, object-oriented interface.

## Workflow Steps

1. **Load data** - Read mzML/mzXML files
2. **Peak detection** - Find chromatographic peaks
3. **RT alignment** - Correct retention time drift
4. **Correspondence** - Group peaks across samples
5. **Gap filling** - Fill missing peak values

## Key Parameters

### CentWave (centroided data)
- `peakwidth`: Expected peak width range (seconds)
- `ppm`: m/z tolerance in ppm
- `snthresh`: Signal-to-noise threshold

### Obiwarp (alignment)
- `binSize`: m/z bin size
- `distFun`: Distance function (cor, cor_opt, cov)

### PeakDensity (grouping)
- `bw`: RT bandwidth for density estimation
- `minFraction`: Min fraction of samples with peak

## Choosing Peak Detection

| Data Type | Method | When to Use |
|-----------|--------|-------------|
| Centroided | CentWave | Most modern instruments |
| Profile | MatchedFilter | Older instruments |
| High-res | CentWave + low ppm | Orbitrap, QTOF |

## QC Recommendations

1. Include pooled QC samples every 10 injections
2. Check TIC for injection issues
3. Verify RT alignment visually
4. PCA should show QC clustering

## References

- XCMS: doi:10.1021/ac051437y
- XCMS3: doi:10.1021/acs.analchem.7b03003
- Documentation: https://bioconductor.org/packages/xcms/
