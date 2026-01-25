---
name: bio-metabolomics-xcms-preprocessing
description: XCMS3 workflow for LC-MS/MS metabolomics preprocessing. Covers peak detection, retention time alignment, correspondence (grouping), and gap filling. Use when processing raw LC-MS data into a feature table for untargeted metabolomics.
tool_type: r
primary_tool: xcms
---

# XCMS Metabolomics Preprocessing

Requires Bioconductor 3.18+ with xcms 4.0+ and MSnbase 2.28+.

## Load Raw Data

```r
library(xcms)
library(MSnbase)

# Read mzML/mzXML files
raw_files <- list.files('raw_data', pattern = '\\.(mzML|mzXML)$', full.names = TRUE)

# Create OnDiskMSnExp object
raw_data <- readMSData(raw_files, mode = 'onDisk')

# Check data
raw_data
table(msLevel(raw_data))
```

## Define Sample Groups

```r
# Sample metadata
sample_info <- data.frame(
    sample_name = basename(raw_files),
    sample_group = c(rep('Control', 5), rep('Treatment', 5), rep('QC', 3)),
    injection_order = 1:length(raw_files)
)

# Assign to phenoData
pData(raw_data) <- sample_info
```

## Peak Detection (Centroided)

```r
# CentWave algorithm for centroided data
cwp <- CentWaveParam(
    peakwidth = c(5, 30),       # Peak width range in seconds
    ppm = 15,                    # m/z tolerance
    snthresh = 10,               # Signal-to-noise threshold
    prefilter = c(3, 1000),      # Min peaks and intensity
    mzdiff = 0.01,               # Minimum m/z difference
    noise = 1000,                # Noise level
    integrate = 1                # Integration method
)

# Run peak detection
xdata <- findChromPeaks(raw_data, param = cwp)

# Summary
head(chromPeaks(xdata))
cat('Peaks found:', nrow(chromPeaks(xdata)), '\n')
```

## Peak Detection (Profile Data)

```r
# MatchedFilter for profile/continuum data
mfp <- MatchedFilterParam(
    binSize = 0.1,
    fwhm = 30,
    snthresh = 10,
    step = 0.1,
    mzdiff = 0.8
)

xdata_profile <- findChromPeaks(raw_data, param = mfp)
```

## Retention Time Alignment

```r
# Obiwarp alignment (recommended)
obp <- ObiwarpParam(
    binSize = 0.5,
    response = 1,
    distFun = 'cor_opt',
    gapInit = 0.3,
    gapExtend = 2.4
)

xdata <- adjustRtime(xdata, param = obp)

# Check alignment
plotAdjustedRtime(xdata)
```

## Peak Correspondence (Grouping)

```r
# Group peaks across samples
pdp <- PeakDensityParam(
    sampleGroups = pData(xdata)$sample_group,
    bw = 5,                      # RT bandwidth
    minFraction = 0.5,           # Min fraction of samples
    minSamples = 1,              # Min samples per group
    binSize = 0.025              # m/z bin size
)

xdata <- groupChromPeaks(xdata, param = pdp)

# Check feature definitions
featureDefinitions(xdata)
cat('Features:', nrow(featureDefinitions(xdata)), '\n')
```

## Gap Filling

```r
# Fill in missing peaks
fpp <- ChromPeakAreaParam()
xdata <- fillChromPeaks(xdata, param = fpp)

# Alternative: FillChromPeaksParam for more control
fpp2 <- FillChromPeaksParam(
    expandMz = 0,
    expandRt = 0,
    ppm = 0
)
```

## Extract Feature Table

```r
# Get feature values (intensity matrix)
feature_values <- featureValues(xdata, method = 'maxint', value = 'into')

# Feature definitions (m/z, RT)
feature_defs <- featureDefinitions(xdata)
feature_defs <- as.data.frame(feature_defs)
feature_defs$feature_id <- rownames(feature_defs)

# Combine
feature_table <- cbind(feature_defs[, c('feature_id', 'mzmed', 'rtmed')], feature_values)
rownames(feature_table) <- feature_table$feature_id

# Save
write.csv(feature_table, 'feature_table.csv', row.names = FALSE)
```

## Quality Control

```r
# TIC for each sample
tic <- chromatogram(raw_data, aggregationFun = 'sum')
plot(tic)

# Peak count per sample
peak_counts <- table(chromPeaks(xdata)[, 'sample'])
barplot(peak_counts, main = 'Peaks per sample')

# Check RT correction
par(mfrow = c(1, 2))
plotAdjustedRtime(xdata, col = pData(xdata)$sample_group)

# PCA of features
library(pcaMethods)
log_values <- log2(feature_values + 1)
log_values[is.na(log_values)] <- 0
pca <- pca(t(log_values), nPcs = 3, method = 'ppca')
plotPcs(pca, col = as.factor(pData(xdata)$sample_group))
```

## CAMERA Annotation (Isotopes/Adducts)

```r
library(CAMERA)

# Create CAMERA object
xsa <- xsAnnotate(as(xdata, 'xcmsSet'))

# Group by RT
xsa <- groupFWHM(xsa, perfwhm = 0.6)

# Find isotopes
xsa <- findIsotopes(xsa, mzabs = 0.01, ppm = 10)

# Find adducts
xsa <- findAdducts(xsa, polarity = 'positive')

# Get annotated peak list
camera_results <- getPeaklist(xsa)
```

## Export for MetaboAnalyst

```r
# Format for MetaboAnalyst web or R package
export_data <- t(feature_values)
colnames(export_data) <- paste0('M', round(feature_defs$mzmed, 4), 'T', round(feature_defs$rtmed, 1))

# Add sample info
export_df <- data.frame(Sample = rownames(export_data), Group = pData(xdata)$sample_group, export_data)

write.csv(export_df, 'metaboanalyst_input.csv', row.names = FALSE)
```

## Related Skills

- metabolite-annotation - Identify metabolites
- normalization-qc - Normalize feature table
- statistical-analysis - Differential analysis
