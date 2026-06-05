#!/usr/bin/env Rscript
# Reference: ichorCNA 0.6.0+, HMMcopy 1.40+ | Verify API if version differs
# ichorCNA tumor-fraction estimation from shallow WGS.
# ichorCNA is a command-line script (scripts/runIchorCNA.R), NOT an importable R function.
# This driver shells out to readCounter (HMMcopy) and runIchorCNA.R, then parses .params.txt.

bin_bam <- function(bam_file, wig_file, window = 1000000, quality = 20,
                    chromosome = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y') {
    cmd <- sprintf('readCounter --window %d --quality %d --chromosome "%s" %s > %s',
                   window, quality, chromosome, bam_file, wig_file)
    system(cmd)
    wig_file
}

run_ichorcna <- function(wig_file, out_dir, ichor_script, gc_wig, map_wig,
                         centromere, normal_panel, sample_id = NULL,
                         normal = 'c(0.5,0.6,0.7,0.8,0.9)', ploidy = 'c(2,3)',
                         max_cn = 7, sc_states = 'c(1,3)', genome_build = 'hg38',
                         genome_style = 'UCSC') {
    if (is.null(sample_id)) {
        sample_id <- gsub('\\.wig$', '', basename(wig_file))
    }
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    args <- c('--id', sample_id, '--WIG', wig_file,
              '--gcWig', gc_wig, '--mapWig', map_wig,
              '--centromere', centromere, '--normalPanel', normal_panel,
              '--normal', shQuote(normal), '--ploidy', shQuote(ploidy),
              '--maxCN', max_cn, '--scStates', shQuote(sc_states),
              '--estimateNormal', 'TRUE', '--estimatePloidy', 'TRUE',  # type=logical: need explicit TRUE
              '--estimateScPrevalence', 'TRUE',
              '--txnE', '0.9999999', '--txnStrength', '1e7',  # segment-length prior (defaults)
              '--minMapScore', '0.9',                         # drop low-mappability bins
              '--genomeBuild', genome_build, '--genomeStyle', genome_style,
              '--outDir', out_dir)
    system2('Rscript', c(ichor_script, args))
    file.path(out_dir, paste0(sample_id, '.params.txt'))
}

# Low-tumor-fraction recipe: seed EM near TF 5/1/0.5/0.1%, fix ploidy=2, drop subclonality.
run_ichorcna_lowtf <- function(wig_file, out_dir, ichor_script, gc_wig, map_wig,
                               centromere, normal_panel, sample_id = NULL,
                               genome_build = 'hg38', genome_style = 'UCSC') {
    run_ichorcna(wig_file, out_dir, ichor_script, gc_wig, map_wig, centromere,
                 normal_panel, sample_id, normal = 'c(0.95,0.99,0.995,0.999)',
                 ploidy = 'c(2)', max_cn = 3, sc_states = 'c()',
                 genome_build = genome_build, genome_style = genome_style)
}

# Tumor fraction = 1 - n_est for the selected (max-loglik) solution, written as row 1.
parse_ichor <- function(params_file) {
    p <- read.table(params_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    data.frame(
        sample = gsub('\\.params\\.txt$', '', basename(params_file)),
        tumor_fraction = 1 - p$n_est[1],
        ploidy = p$phi_est[1],
        loglik = p$loglik[1]
    )
}

parse_cohort <- function(results_dir) {
    files <- list.files(results_dir, pattern = '\\.params\\.txt$', full.names = TRUE, recursive = TRUE)
    do.call(rbind, lapply(files, parse_ichor))
}

if (sys.nframe() == 0) {
    cat('ichorCNA tumor-fraction estimation\n')
    cat('bin_bam() -> run_ichorcna() [or run_ichorcna_lowtf()] -> parse_ichor()\n')
    cat('parse_cohort() summarises a results directory; TF = 1 - n_est\n')
}
