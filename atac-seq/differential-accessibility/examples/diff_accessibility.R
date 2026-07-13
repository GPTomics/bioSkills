#!/usr/bin/env Rscript
# Reference: DiffBind 3.12+, DESeq2 1.42+, edgeR 4.0+, ChIPseeker 1.38+, sva 3.50+ | Verify API if version differs
# DiffBind workflow with summits=250 fixed-width counting and ENCODE-style consensus.
# Includes optional SVA hidden-batch correction and library-size normalization for global-shift biology.

suppressPackageStartupMessages({
    library(DiffBind); library(rtracklayer); library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(sva)
})

run_diff <- function(sample_sheet, output_prefix = 'diff_atac',
                     normalize_mode = DBA_NORM_NATIVE,    # Switch to DBA_NORM_LIB for global shift
                     fdr_thr = 0.05, lfc_thr = 1,         # Standard ENCODE-style thresholds
                     min_overlap = 2,                     # Majority rule: peak in >= 2 reps
                     fixed_width_summit = 250,            # +/- 250 = 501 bp consensus (Corces 2018)
                     use_sva = FALSE) {
    cat('=== DiffBind workflow ===\n')
    dba <- dba(sampleSheet = sample_sheet)

    # Consensus peakset: re-center on summits +/- 250 bp; require peak in >= min_overlap reps
    cat(sprintf('Counting reads in fixed-width consensus (summits +/- %d, minOverlap=%d)...\n',
                fixed_width_summit, min_overlap))
    dba <- dba.count(dba, summits = fixed_width_summit, minOverlap = min_overlap, bParallel = TRUE)

    # Normalize: NATIVE = per-tool native method (RLE for DESeq2, TMM for edgeR); LIB = full library size for global shifts.
    # DiffBind's own default is DBA_NORM_LIB; reads-in-peaks is a library-size option (library = DBA_LIBSIZE_PEAKREADS).
    cat('Normalizing (', deparse(substitute(normalize_mode)), ')...\n')
    dba <- dba.normalize(dba, normalize = normalize_mode)

    # Optional: hidden-batch SVA before fitting
    if (use_sva) {
        cat('Estimating surrogate variables with SVA...\n')
        se <- dba(dba, bSummarizedExperiment = TRUE)
        counts_mat <- SummarizedExperiment::assay(se)
        meta <- as.data.frame(dba(dba)$samples)
        mod  <- model.matrix(~Condition, meta)
        mod0 <- model.matrix(~1, meta)
        sv <- svaseq(counts_mat[rowMeans(counts_mat) > 1, ], mod, mod0, n.sv = 2)$sv
        # SVA covariates would feed back via dba.contrast(..., design=...) custom -- example shows extraction
        cat(sprintf('  -> SV1, SV2 estimated (rank=%d)\n', ncol(sv)))
    }

    cat('Setting up contrast...\n')
    dba <- dba.contrast(dba, contrast = c('Condition', 'treated', 'control'),
                        minMembers = 2)
    cat('Fitting DESeq2 model...\n')
    dba <- dba.analyze(dba, method = DBA_DESEQ2)

    # Report applies threshold within DiffBind; DBA reporting columns: Conc, Conc_X, Conc_Y, Fold, p.value, FDR
    cat(sprintf('Reporting at FDR < %.2f, |Fold| >= %.1f\n', fdr_thr, lfc_thr))
    results <- dba.report(dba, th = fdr_thr, fold = lfc_thr)
    cat(sprintf('  Differentially accessible peaks: %d\n', length(results)))
    cat(sprintf('  Opened (Fold > 0): %d   Closed (Fold < 0): %d\n',
                sum(results$Fold > 0), sum(results$Fold < 0)))

    # Export BEDs (gain / loss)
    export(results[results$Fold > 0], sprintf('%s_opened.bed', output_prefix))
    export(results[results$Fold < 0], sprintf('%s_closed.bed', output_prefix))

    # Annotate to genes (promoter region -2kb to +500bp; ChIPseeker default is too wide)
    cat('Annotating peaks to genes...\n')
    peakAnno <- annotatePeak(results, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                             tssRegion = c(-2000, 500), level = 'gene', verbose = FALSE)
    write.csv(as.data.frame(peakAnno), sprintf('%s_annotated.csv', output_prefix), row.names = FALSE)
    pdf(sprintf('%s_annoplot.pdf', output_prefix))
    plotAnnoPie(peakAnno); plotDistToTSS(peakAnno)
    dev.off()

    # Diagnostic plots
    pdf(sprintf('%s_diagnostics.pdf', output_prefix), width = 8, height = 6)
    dba.plotPCA(dba, attributes = DBA_CONDITION, label = DBA_ID)
    dba.plotMA(dba); dba.plotVolcano(dba); dba.plotHeatmap(dba, contrast = 1, correlations = FALSE)
    dev.off()

    invisible(list(dba = dba, results = results, anno = peakAnno))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) run_diff(args[1])
