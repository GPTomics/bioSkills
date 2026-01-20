#!/usr/bin/env Rscript
# Differential accessibility with DiffBind

library(DiffBind)

run_diff_accessibility <- function(sample_sheet, output_prefix = 'diff_atac') {
    cat('Loading samples...\n')
    dba <- dba(sampleSheet = sample_sheet)

    cat('Counting reads in peaks...\n')
    dba <- dba.count(dba)

    cat('Setting up contrast...\n')
    dba <- dba.contrast(dba, categories = DBA_CONDITION)

    cat('Running differential analysis...\n')
    dba <- dba.analyze(dba)

    cat('Extracting results...\n')
    results <- dba.report(dba)

    # Export results
    export(results, paste0(output_prefix, '_results.bed'))

    # Summary
    cat(sprintf('\nDifferentially accessible regions: %d\n', length(results)))
    cat(sprintf('More accessible: %d\n', sum(results$Fold > 0)))
    cat(sprintf('Less accessible: %d\n', sum(results$Fold < 0)))

    # MA plot
    pdf(paste0(output_prefix, '_MA.pdf'))
    dba.plotMA(dba)
    dev.off()

    return(results)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
    run_diff_accessibility(args[1])
}
