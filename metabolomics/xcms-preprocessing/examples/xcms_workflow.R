# Reference: xcms 4.x+, MsExperiment 1.4+, faahKO 1.42+ | Verify API if version differs
# End-to-end modern-API xcms preprocessing on the bundled faahKO test data:
# raw -> peak detection -> alignment -> correspondence -> gap-filling -> feature table.
# Output is written to tempdir() and removed; no artifacts are left in the repo.

library(xcms)
library(MsExperiment)
library(faahKO)

cdfs <- dir(system.file('cdf', package = 'faahKO'), full.names = TRUE, recursive = TRUE)[c(1, 2, 7, 8)]
pd <- data.frame(sample_name = sub('.CDF', '', basename(cdfs), fixed = TRUE),
                 sample_group = c('KO', 'KO', 'WT', 'WT'))

raw <- readMsExperiment(spectraFiles = cdfs, sampleData = pd)
cat('Loaded', length(cdfs), 'files\n')

# CentWave for centroided high-res data. faahKO is broad-peak conventional-HPLC LC-MS data,
# so peakwidth runs broad (20-80 s); for UHPLC measure EIC base-widths and narrow it.
# ppm reflects across-scan centroid scatter, not the spec mass accuracy.
cwp <- CentWaveParam(ppm = 25, peakwidth = c(20, 80), snthresh = 10,
                     prefilter = c(3, 1000), noise = 1000)
xdata <- findChromPeaks(raw, param = cwp)
cat('Peaks detected:', nrow(chromPeaks(xdata)), '\n')

# obiwarp needs no prior peaks; binSize here is the m/z profile bin (distinct from grouping binSize).
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))

# Peak-density correspondence. bw should track residual post-alignment RT scatter;
# minFraction 0.5 keeps features present in at least half of one group.
pdp <- PeakDensityParam(sampleGroups = sampleData(xdata)$sample_group,
                        bw = 30, minFraction = 0.5, binSize = 0.025)
xdata <- groupChromPeaks(xdata, param = pdp)
cat('Features:', nrow(featureDefinitions(xdata)), '\n')

# Gap-filling fabricates intensities for absent compounds, so track the filled flag.
xdata <- fillChromPeaks(xdata, param = ChromPeakAreaParam())
filled_fraction <- mean(chromPeakData(xdata)$is_filled)
cat('Fraction of filled peaks:', round(filled_fraction, 3), '\n')

feat <- featureValues(xdata, value = 'into')
defs <- as.data.frame(featureDefinitions(xdata))
result <- data.frame(feature = rownames(feat), mz = defs$mzmed, rt = defs$rtmed, feat,
                     row.names = NULL)

out <- file.path(tempdir(), 'feature_table.csv')
write.csv(result, out, row.names = FALSE)
cat('Wrote', nrow(result), 'features to', out, '\n')
file.remove(out)
