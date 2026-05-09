#!/usr/bin/env Rscript
# Reference: FRASER 2.0 (>=1.99.0), R 4.4+, BiocManager 1.30+ | Verify API if version differs
# FRASER 2.0 outlier splicing detection for rare-disease diagnostics
#
# Workflow:
# 1. Build FraserDataSet from patient + control BAMs
# 2. Count split reads per junction
# 3. Compute Intron Jaccard Index (FRASER 2.0 default metric)
# 4. Fit Beta-binomial autoencoder; flag outliers
# 5. Filter aberrant junctions in patient sample

suppressPackageStartupMessages({
    library(FRASER)
    library(BiocParallel)
    library(dplyr)
})

bam_dir <- 'bams'
working_dir <- 'fraser_workdir'
patient_id <- 'PATIENT_001'

bam_files <- list.files(bam_dir, pattern = '\\.bam$', full.names = TRUE)
sample_table <- data.frame(
    sampleID = gsub('\\.bam$', '', basename(bam_files)),
    bamFile = bam_files,
    pairedEnd = TRUE
)

settings <- FraserDataSet(
    colData = sample_table,
    workingDir = working_dir,
    name = 'rare_disease_cohort'
)

settings <- countRNAData(settings, BPPARAM = MulticoreParam(8))

fds <- calculatePSIValues(settings)

fds <- filterExpressionAndVariability(
    fds,
    minDeltaPsi = 0.0,
    minExpressionInOneSample = 20,
    quantile = 0.05,
    quantileMinExpression = 1
)

# FRASER 2.0: use Intron Jaccard Index
fitMetrics(fds) <- 'jaccard'

# Tune q hyperparameter once per cohort; q=10 default
fds <- FRASER(
    fds,
    q = c(jaccard = 10),
    BPPARAM = MulticoreParam(8)
)

# Default delta cutoff in FRASER 2.0 is 0.1 (was 0.3 in v1.x)
all_results <- results(
    fds,
    psiType = 'jaccard',
    padjCutoff = 0.05,
    deltaPsiCutoff = 0.1
)

patient_results <- all_results %>%
    as.data.frame() %>%
    filter(sampleID == patient_id) %>%
    arrange(padjust)

write.table(
    patient_results,
    file = file.path(working_dir, paste0(patient_id, '_outliers.tsv')),
    sep = '\t', quote = FALSE, row.names = FALSE
)

# Plot diagnostic for this patient
plotVolcano(fds, sampleID = patient_id, type = 'jaccard')

cat(sprintf('%d aberrant junctions in %s (padj<0.05, |delta|>=0.1)\n',
            nrow(patient_results), patient_id))
