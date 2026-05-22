#!/usr/bin/env Rscript
# Reference: scarHRD 0.1.1+ | Verify API if version differs
#
# Compute the three HRD genomic scars (LOH, LST, TAI) and their sum from allele-specific
# copy number. scarHRD needs allele-specific input -- a Sequenza .seqz file or an
# ASCAT-style A/B-allele segment table -- NOT relative log2 copy ratio.
#
# LST rises with ploidy: a whole-genome-doubled tumor scores falsely high unless the
# score is computed with the correct ploidy. This script reports the components so an
# LST-driven high score can be cross-checked against WGD status.

suppressMessages(library(scarHRD))

args <- commandArgs(trailingOnly = TRUE)
infile    <- ifelse(length(args) >= 1, args[1], 'sample.small.seqz.gz')
reference <- ifelse(length(args) >= 2, args[2], 'grch38')
# is_seqz=TRUE for a Sequenza .seqz file; FALSE for a pre-computed allele-specific table.
is_seqz   <- ifelse(length(args) >= 3, as.logical(args[3]), TRUE)

hrd <- scar_score(infile, reference = reference, seqz = is_seqz)

# hrd columns (note the space in 'Telomeric AI' and hyphen in 'HRD-sum'):
# 'HRD' (LOH count), 'Telomeric AI', 'LST', 'HRD-sum'.
print(hrd)

loh <- hrd[['HRD']]
tai <- hrd[['Telomeric AI']]
lst <- hrd[['LST']]
total <- hrd[['HRD-sum']]
cat(sprintf('\nLOH=%d  TAI=%d  LST=%d  HRD-sum=%d\n', loh, tai, lst, total))

# If LST dominates the sum, suspect whole-genome doubling rather than true HR deficiency.
if (lst > (loh + tai))
    cat('NOTE: LST dominates the score. Confirm the ploidy used is correct and check\n',
        '      whole-genome-doubling status before calling this tumor HRD-positive.\n')

cat('\nReminder: the HRD score is a SCAR of PAST HR deficiency. A BRCA-reverted tumor\n',
    'still scores high. Integrate current HR-pathway status before predicting response.\n')
cat('Use the cutoff validated for the specific assay (e.g. GIS >= 42 for Myriad\n',
    'myChoice); cutoffs are not portable across assays.\n')
