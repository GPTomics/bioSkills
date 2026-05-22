#!/usr/bin/env Rscript
# Reference: FACETS 0.6+ (R package + snp-pileup) | Verify API if version differs
#
# Two-pass FACETS allele-specific copy number for a tumor-normal pair.
# The two-pass design (coarse purity run -> dipLogR-seeded sensitivity run) stabilizes
# the fit against the purity-ploidy identifiability problem. snp-pileup must be run
# first to produce the input CSV:
#   snp-pileup -g -q15 -Q20 -P100 -r25,0 dbsnp_common.vcf.gz \
#       sample.snp_pileup.csv.gz normal.bam tumor.bam

suppressMessages(library(facets))

args <- commandArgs(trailingOnly = TRUE)
pileup    <- ifelse(length(args) >= 1, args[1], 'sample.snp_pileup.csv.gz')
sample_id <- ifelse(length(args) >= 2, args[2], 'sample')
gbuild    <- ifelse(length(args) >= 3, args[3], 'hg38')

# cval: segmentation critical value. Panels/WES ~150-300; WGS ~25-100.
# Too low causes hyperfragmentation (spurious micro-segments).
cval_purity   <- 300   # coarse pass: stable purity/ploidy
cval_focal    <- 150   # fine pass: focal sensitivity

set.seed(1234)         # FACETS uses random initialization; fix the seed for reproducibility
rcmat <- readSnpMatrix(pileup)
xx <- preProcSample(rcmat, gbuild = gbuild)

# Pass 1: estimate purity, ploidy, and the diploid baseline (dipLogR).
oo1 <- procSample(xx, cval = cval_purity)
fit1 <- emcncf(oo1)

# Pass 2: re-segment at finer cval, seeded by pass-1 dipLogR so focal events are
# recovered without re-litigating the diploid baseline.
oo2 <- procSample(xx, cval = cval_focal, dipLogR = oo1$dipLogR)
fit2 <- emcncf(oo2)

cat(sprintf('%s: purity=%.3f ploidy=%.3f dipLogR=%.3f\n',
            sample_id, fit2$purity, fit2$ploidy, oo2$dipLogR))

# cncf columns: tcn.em = total copy number, lcn.em = minor copy number, cf.em = cell
# fraction. lcn.em == 0 marks loss of heterozygosity.
cncf <- fit2$cncf
cncf$loh <- !is.na(cncf$lcn.em) & cncf$lcn.em == 0
write.table(cncf, file = paste0(sample_id, '.facets.cncf.tsv'),
            sep = '\t', quote = FALSE, row.names = FALSE)

# ALWAYS inspect the diagnostic plot. A jagged profile or tcn.em incoherent with
# cnlr.median means the fit is bad -- re-tune cval.
pdf(paste0(sample_id, '.facets.pdf'), width = 9, height = 7)
plotSample(x = oo2, emfit = fit2)
dev.off()

if (fit2$purity < 0.40)
    cat('WARNING: purity < 0.40 -- below the reliable range for allele-specific calling.\n')
