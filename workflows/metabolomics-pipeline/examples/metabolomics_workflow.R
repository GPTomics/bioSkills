# Reference: xcms 4.x+, pmp 1.14+, ropls 1.34+, MetaboAnalystR 4.0+ | Verify API if version differs
#
# End-to-end untargeted LC-MS metabolomics orchestration. Stage 1 (raw mzML -> feature
# table) is shown with the modern xcms 4.x MsExperiment API but commented out, because it
# needs real centroided mzML files. Stages 2-5 run on a synthetic feature table so the QC /
# statistics / pathway logic is demonstrable without instrument data. Each stage defers to
# its component skill for parameter rationale and failure modes.

# === STAGE 1: FEATURE EXTRACTION (xcms 4.x) -- see metabolomics/xcms-preprocessing ===
# pd is a data.frame, one row per file, with a sample_group column.
#
# library(xcms)
# raw <- readMsExperiment(spectraFiles = mzml_files, sampleData = pd)
# cwp <- CentWaveParam(ppm = 10, peakwidth = c(2, 20), snthresh = 10,
#                      prefilter = c(3, 1000), noise = 1000)
# xdata <- findChromPeaks(raw, param = cwp)
# xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))
# pdp <- PeakDensityParam(sampleGroups = sampleData(xdata)$sample_group,
#                         bw = 5, minFraction = 0.5, binSize = 0.025)
# xdata <- groupChromPeaks(xdata, param = pdp)          # group on corrected RT (obiwarp needs no pre-grouping)
# xdata <- fillChromPeaks(xdata, param = ChromPeakAreaParam())
# feat <- featureValues(xdata, value = 'into')          # features x samples; filled cells are imputations
# defs <- featureDefinitions(xdata)                     # mzmed / rtmed per feature

output_dir <- file.path(tempdir(), 'metabolomics_demo')
dir.create(output_dir, showWarnings = FALSE)

# === SYNTHETIC FEATURE TABLE standing in for Stage 1 output ===
set.seed(1)
n_features <- 300
n_per_group <- 12
n_qc <- 6
sample_names <- c(paste0('QC', 1:n_qc), paste0('Ctrl', 1:n_per_group), paste0('Trt', 1:n_per_group))
sample_group <- c(rep('QC', n_qc), rep('Control', n_per_group), rep('Treatment', n_per_group))
batch_id <- rep(1, length(sample_names))

base_intensity <- 2^runif(n_features, 8, 20)
# QCs are matrix-matched pools: tight technical noise. Study samples carry biological spread,
# so a real feature has small QC SD relative to biological SD (low D-ratio) by construction.
qc_noise <- matrix(rlnorm(n_features * n_qc, 0, 0.05), nrow = n_features)
bio_noise <- matrix(rlnorm(n_features * (2 * n_per_group), 0, 0.20), nrow = n_features)
feat <- base_intensity * cbind(qc_noise, bio_noise)

# Inject a true Treatment effect into 20 features (3-fold up -> log2FC ~1.6), clearly ABOVE the
# |log2FC|>1 significance cutoff so recovery measures sensitivity, not cutoff-boundary collision.
true_hits <- 1:20
feat[true_hits, sample_group == 'Treatment'] <- feat[true_hits, sample_group == 'Treatment'] * 3

# BLOCK-RANDOMIZE run order: study samples get interleaved random run positions so group is NOT
# collinear with injection order/drift -- the "original sin" (commitment #3) an ordered layout bakes in.
# (Defined here, after the noise draws, so it does not perturb the RNG stream above.)
study_idx <- which(sample_group != 'QC')
injection_order <- integer(length(sample_names))
injection_order[sample_group == 'QC'] <- seq_len(n_qc)                       # QCs condition the column first
injection_order[study_idx] <- n_qc + sample(seq_along(study_idx))            # then Control/Treatment in random order

# Inject injection-order drift (uniform per-sample sensitivity loss); PQN in Stage 2 removes this
# uniform per-column scaling -- but only because order is randomized, not confounded with group.
drift <- 1 - 0.3 * (injection_order / max(injection_order))
feat <- sweep(feat, 2, drift, '*')

rownames(feat) <- paste0('FT', sprintf('%03d', 1:n_features))
colnames(feat) <- sample_names
mz <- runif(n_features, 80, 900)

# === STAGE 2: QC, DRIFT, NORMALIZATION -- see metabolomics/normalization-qc ===
# Real pipeline: filter_peaks_by_fraction -> QCRSC -> filter_peaks_by_rsd -> pqn_normalisation
# (pmp; features in ROWS). Demonstrated here with base R so the script runs without Bioconductor.

is_qc <- sample_group == 'QC'

qc_rsd <- apply(feat[, is_qc], 1, function(x) sd(x) / mean(x))
sd_qc <- apply(feat[, is_qc], 1, sd)
sd_bio <- apply(feat[, !is_qc], 1, sd)
dratio <- sd_qc / sd_bio
keep <- qc_rsd <= 0.30 & dratio <= 0.50            # Broadhurst 2018: RSD<=30%, D-ratio<=0.5
feat_qc <- feat[keep, ]
cat('Features kept after RSD/D-ratio filter:', nrow(feat_qc), 'of', n_features, '\n')

# PQN: median-quotient against the QC reference spectrum, robust to a minority of changed features.
reference <- apply(feat_qc[, is_qc], 1, median)
quotients <- sweep(feat_qc, 1, reference, '/')
dilution_factor <- apply(quotients, 2, median, na.rm = TRUE)
normalized <- sweep(feat_qc, 2, dilution_factor, '/')
logged <- log2(normalized)

# === STAGE 4: STATISTICS -- see metabolomics/statistical-analysis ===
# Univariate Welch + BH FDR on the study samples (QCs excluded). A permutation-validated
# OPLS-DA via ropls::opls(permI=1000) is the multivariate half; shown in SKILL.md.
study <- sample_group %in% c('Control', 'Treatment')
study_group <- factor(sample_group[study])

welch <- apply(logged[, study], 1, function(x)
    t.test(x[study_group == 'Treatment'], x[study_group == 'Control'])$p.value)
lfc <- apply(logged[, study], 1, function(x)
    mean(x[study_group == 'Treatment']) - mean(x[study_group == 'Control']))
padj <- p.adjust(welch, method = 'BH')             # default p.adjust is 'holm'; BH is the discovery default

results <- data.frame(feature_id = rownames(logged), mz = mz[keep],
                      log2fc = lfc, pval = welch, padj = padj)
results$significant <- results$padj < 0.05 & abs(results$log2fc) > 1   # FDR 5% + 2-fold
cat('Significant features (FDR<0.05, |log2FC|>1):', sum(results$significant), '\n')
recovered <- sum(results$significant & results$feature_id %in% rownames(feat)[true_hits])
cat('True planted hits recovered:', recovered, 'of', length(true_hits), '\n')

# === STAGE 3 + 5: ANNOTATION & PATHWAY MAPPING -- see metabolite-annotation, pathway-mapping ===
# Stage 3 attaches an MSI/Schymanski confidence level to each significant feature (no level,
# no identification). Stage 5 then chooses ORA (identified compounds, assay-coverage background)
# or mummichog (raw m/z, FULL feature table as background). Both deferred to their skills:
#
# library(MetaboAnalystR)
# mSet <- InitDataObjects('mass_all', 'mummichog', FALSE)
# mSet <- UpdateInstrumentParameters(mSet, 5.0, 'negative')   # ppm + ionization mode mandatory
# mSet <- Read.PeakListData(mSet, 'peaks_full.txt')           # ENTIRE table, not significant-only
# mSet <- PerformPSEA(mSet, 'hsa_mfn', 'current', permNum = 1000)

# === OUTPUT to tempdir, then clean up ===
results_path <- file.path(output_dir, 'differential_metabolites.csv')
matrix_path <- file.path(output_dir, 'normalized_feature_matrix.csv')
write.csv(results, results_path, row.names = FALSE)
write.csv(normalized, matrix_path)
cat('Wrote results to', output_dir, '\n')

unlink(output_dir, recursive = TRUE)
cat('Cleaned up temporary outputs\n')
