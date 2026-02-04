#!/usr/bin/env Rscript
# Tumor fraction estimation with ichorCNA

library(ichorCNA)

run_ichorcna <- function(wig_file, output_dir, sample_id = NULL,
                         gc_wig, map_wig, normal_panel, centromere) {
    if (is.null(sample_id)) {
        sample_id <- gsub('.wig$', '', basename(wig_file))
    }

    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    runIchorCNA(
        WIG = wig_file,
        gcWig = gc_wig,
        mapWig = map_wig,
        normalPanel = normal_panel,
        centromere = centromere,
        outDir = output_dir,
        id = sample_id,
        normal = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99),
        ploidy = c(2, 3),
        maxCN = 5,
        estimateScPrevalence = TRUE,
        scStates = c(1, 3),
        txnE = 0.9999,
        txnStrength = 10000,
        chrs = paste0('chr', c(1:22, 'X'))
    )

    return(file.path(output_dir, paste0(sample_id, '.params.txt')))
}

parse_results <- function(results_dir) {
    param_files <- list.files(results_dir, pattern = '.params.txt$',
                              full.names = TRUE, recursive = TRUE)

    results <- data.frame()

    for (f in param_files) {
        params <- read.table(f, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
        sample_id <- gsub('.params.txt$', '', basename(f))

        results <- rbind(results, data.frame(
            sample = sample_id,
            tumor_fraction = 1 - params$n[1],
            ploidy = params$phi[1],
            log_likelihood = params$loglik[1]
        ))
    }

    return(results)
}

batch_process <- function(wig_dir, output_dir, gc_wig, map_wig, normal_panel, centromere, n_cores = 4) {
    library(parallel)

    wig_files <- list.files(wig_dir, pattern = '.wig$', full.names = TRUE)

    process_one <- function(wig_file) {
        tryCatch({
            run_ichorcna(wig_file, output_dir, gc_wig = gc_wig,
                        map_wig = map_wig, normal_panel = normal_panel,
                        centromere = centromere)
            return(list(file = wig_file, status = 'success'))
        }, error = function(e) {
            return(list(file = wig_file, status = 'failed', error = e$message))
        })
    }

    results <- mclapply(wig_files, process_one, mc.cores = n_cores)
    return(results)
}

cat('ichorCNA Tumor Fraction Estimation\n')
cat('===================================\n\n')
cat('Usage:\n')
cat('1. run_ichorcna() - Process single sample\n')
cat('2. batch_process() - Process multiple samples\n')
cat('3. parse_results() - Extract tumor fractions\n')
