#!/usr/bin/env Rscript
# Reference: DNAcopy 1.76+ | Verify API if version differs
#
# Circular binary segmentation of a normalized log2 copy-ratio profile.
# Segments per chromosome (CBS otherwise bridges assembly gaps) and merges adjacent
# segments within a noise-scaled threshold (undo.SD) to guard against oversegmentation.
# Input: a TSV with columns chrom, maploc (bin midpoint), log2.

suppressMessages(library(DNAcopy))

args <- commandArgs(trailingOnly = TRUE)
infile    <- ifelse(length(args) >= 1, args[1], 'bins.tsv')
sample_id <- ifelse(length(args) >= 2, args[2], 'sample')
outfile   <- ifelse(length(args) >= 3, args[3], paste0(sample_id, '.segments.tsv'))

# alpha:    breakpoint significance (lower = fewer, more confident breakpoints)
# undo.SD:  merge segments whose means are within this many noise SD of each other --
#           the main control against oversegmentation
ALPHA   <- 0.01
UNDO_SD <- 2

bins <- read.delim(infile)
stopifnot(all(c('chrom', 'maploc', 'log2') %in% colnames(bins)))

cna <- CNA(genomdat = bins$log2, chrom = bins$chrom, maploc = bins$maploc,
           data.type = 'logratio', sampleid = sample_id)

# smooth.CNA damps single-bin outliers so one spiky bin does not force a breakpoint.
cna <- smooth.CNA(cna)

seg <- segment(cna, alpha = ALPHA, undo.splits = 'sdundo', undo.SD = UNDO_SD,
               verbose = 1)

out <- seg$output
write.table(out, outfile, sep = '\t', quote = FALSE, row.names = FALSE)

# Oversegmentation check: if segment count is implausibly high relative to genome size,
# the input is too noisy or alpha/undo.SD are too liberal.
n_seg <- nrow(out)
cat(sprintf('%s: %d segments\n', sample_id, n_seg))
if (n_seg > 2000)
    cat('WARNING: high segment count -- suspect oversegmentation. Raise undo.SD or',
        'denoise the input before trusting these segments.\n')

# Reminder: CBS gives RELATIVE copy ratio. Before any gain/loss interpretation, confirm
# the diploid baseline is correctly anchored (see allele-specific-copy-number for an
# absolute ploidy estimate). Centering on the data mode inverts calls in WGD genomes.
