library(xcms)
library(MSnbase)

# Configuration
raw_dir <- 'raw_data'
output_file <- 'feature_table.csv'

# Load data
raw_files <- list.files(raw_dir, pattern = '\\.mzML$', full.names = TRUE)
raw_data <- readMSData(raw_files, mode = 'onDisk')

# Sample info (modify based on your experiment)
pData(raw_data)$sample_group <- ifelse(grepl('ctrl', basename(raw_files)), 'Control', 'Treatment')

cat('Loaded', length(raw_files), 'files\n')

# Peak detection (CentWave)
# peakwidth: Expected chromatographic peak width in seconds; typical 5-30 for LC-MS
# ppm: Mass accuracy; use 5-15 for Orbitrap, 20-30 for TOF
# snthresh: Signal-to-noise; 10 is standard, lower for weak signals
cwp <- CentWaveParam(peakwidth = c(5, 30), ppm = 15, snthresh = 10, prefilter = c(3, 1000))
xdata <- findChromPeaks(raw_data, param = cwp)
cat('Peaks detected:', nrow(chromPeaks(xdata)), '\n')

# RT alignment (Obiwarp)
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.5))
cat('RT aligned\n')

# Correspondence (grouping)
pdp <- PeakDensityParam(sampleGroups = pData(xdata)$sample_group, bw = 5, minFraction = 0.5)
xdata <- groupChromPeaks(xdata, param = pdp)
cat('Features:', nrow(featureDefinitions(xdata)), '\n')

# Gap filling
xdata <- fillChromPeaks(xdata, param = ChromPeakAreaParam())

# Extract feature table
features <- featureValues(xdata, method = 'maxint', value = 'into')
feature_info <- as.data.frame(featureDefinitions(xdata))

result <- data.frame(
    feature = rownames(features),
    mz = feature_info$mzmed,
    rt = feature_info$rtmed,
    features
)

write.csv(result, output_file, row.names = FALSE)
cat('Saved to', output_file, '\n')
