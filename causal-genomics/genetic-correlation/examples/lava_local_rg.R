# Reference: LDSC v1.0.1+, HDL, LAVA | Verify API if version differs
# LAVA local genetic correlation per LDetect locus.
# Detects shared etiology hidden by global rg cancellation.
#
# Operational rule (Werme 2022 Nat Genet 54:274): filter loci on
# univariate local h2 p < 0.05/N_loci in BOTH traits BEFORE running
# bivariate rg. Skipping the filter produces spurious +/- 1 boundary
# estimates at unidentified loci.

suppressPackageStartupMessages({
    library(LAVA)
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
info_file <- if (length(args) >= 1) args[1] else 'input.info.txt'
overlap_file <- if (length(args) >= 2) args[2] else 'sample.overlap.txt'
ref_prefix <- if (length(args) >= 3) args[3] else '1kg_EUR_chr'
loci_file <- if (length(args) >= 4) args[4] else 'blocks_s2500_m25_f1_w200.GRCh38.locfile'
phenos <- if (length(args) >= 5) strsplit(args[5], ',')[[1]] else c('trait1', 'trait2')
out_prefix <- if (length(args) >= 6) args[6] else 'lava_local_rg'

input <- process.input(
    input.info.file = info_file,
    sample.overlap.file = overlap_file,
    ref.prefix = ref_prefix,
    phenos = phenos
)

loci <- read.loci(loci_file)
N_loci <- nrow(loci)
message('LDetect loci: ', N_loci)

univ_alpha <- 0.05 / N_loci  # Bonferroni for univariate local h2 (per Werme 2022)
biv_alpha <- 0.05 / N_loci  # Bonferroni for bivariate local rg

univ_rows <- list()
biv_rows <- list()

for (i in seq_len(N_loci)) {
    locus <- tryCatch(process.locus(loci[i, ], input), error = function(e) NULL)
    if (is.null(locus)) next
    univ <- run.univ(locus)
    univ$locus_id <- loci$LOC[i]
    univ_rows[[length(univ_rows) + 1L]] <- univ

    pass_univ <- all(univ$p < univ_alpha)
    if (!pass_univ) next

    biv <- tryCatch(run.bivar(locus), error = function(e) NULL)
    if (is.null(biv)) next
    biv$locus_id <- loci$LOC[i]
    biv_rows[[length(biv_rows) + 1L]] <- biv
}

univ_df <- rbindlist(univ_rows, fill = TRUE)
biv_df <- rbindlist(biv_rows, fill = TRUE)

if (nrow(biv_df) > 0) {
    biv_df[, padj := p.adjust(p, method = 'bonferroni', n = N_loci)]
    sig_loci <- biv_df[padj < 0.05]
    message('Bonferroni-significant local rg loci: ', nrow(sig_loci))
} else {
    sig_loci <- biv_df
    message('No loci passed univariate h2 filter in both traits.')
}

fwrite(univ_df, paste0(out_prefix, '.univ.tsv'), sep = '\t')
fwrite(biv_df, paste0(out_prefix, '.biv.tsv'), sep = '\t')
fwrite(sig_loci, paste0(out_prefix, '.sig.tsv'), sep = '\t')

if (nrow(sig_loci) > 0) {
    message('Top local-rg hits:')
    print(sig_loci[order(p)][1:min(10, .N)])
}
