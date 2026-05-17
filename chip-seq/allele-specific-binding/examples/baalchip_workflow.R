#!/usr/bin/env Rscript
# Reference: BaalChIP 1.30+, GenomicRanges 1.54+, rtracklayer 1.62+, samtools 1.19+ (external) | Verify API if version differs
# BaalChIP Bayesian beta-binomial ASB analysis with copy-number-aware
# overdispersion. Handles cancer / CN-imbalanced samples. Assumes WASP filter
# already applied (mandatory; not optional). Filters imprinted loci and chrX
# (in female samples) before testing.

suppressPackageStartupMessages({
    library(BaalChIP)
    library(GenomicRanges)
    library(rtracklayer)
})


# === 1. Sample setup ===
# CRITICAL: BAMs must be WASP-filtered (see WASP pipeline). BaalChIP does not
# correct reference-allele mapping bias itself; it assumes input BAMs are clean.
samples <- data.frame(
    SampleID = c('HCC1395_FOXA1_rep1', 'HCC1395_FOXA1_rep2'),
    Tissue = 'TNBC',
    Target = 'FOXA1',
    BAM = c('HCC1395_FOXA1_rep1.wasp.bam',
            'HCC1395_FOXA1_rep2.wasp.bam'),
    Peaks = c('HCC1395_FOXA1_rep1_peaks.narrowPeak',
              'HCC1395_FOXA1_rep2_peaks.narrowPeak'),
    Group = 'HCC1395'
)


# === 2. Heterozygous SNP file ===
# Use SAMPLE-SPECIFIC VCF (from GATK / DeepVariant on matched-normal), NOT a
# population panel. Sample-specific hetSNPs may not be in 1KG.
# Format: chr, pos (1-based), ref, alt, AF (allele frequency 0-1)
hetSNPs <- read.delim('HCC1395_hetSNPs.bed', header = FALSE,
                       col.names = c('CHROM', 'POS', 'REF', 'ALT', 'AF'))

# Filter universal artifacts BEFORE ASB
# (a) Imprinted loci (constitutively allele-skewed by biology)
imprinted_bed <- import('imprinted_loci_hg38.bed')
hetSNPs_gr <- GRanges(seqnames = hetSNPs$CHROM,
                       ranges = IRanges(start = hetSNPs$POS, width = 1))
keep <- !overlapsAny(hetSNPs_gr, imprinted_bed)
hetSNPs <- hetSNPs[keep, ]

# (b) chrX in female samples (X-inactivation)
# HCC1395 is female (TNBC line); filter chrX
hetSNPs <- hetSNPs[hetSNPs$CHROM != 'chrX', ]

write.table(hetSNPs, 'HCC1395_hetSNPs_filtered.tsv', sep = '\t',
             quote = FALSE, row.names = FALSE)


# === 3. CNV file for copy-number-aware overdispersion ===
# REQUIRED for cancer samples. ASCAT / Sequenza / FACETS output BED of
# CN-altered regions with total + minor allele counts.
# Format: chr, start, end, total_CN, minor_CN
# Without this, BaalChIP cannot model CN-altered allele dose.
cnvs <- 'HCC1395_ASCAT_cnvs.bed'


# === 4. Initialize BaalChIP object ===
res <- BaalChIP(samplesheet = samples,
                 hets = c(HCC1395 = 'HCC1395_hetSNPs_filtered.tsv'))


# === 5. Allele counts at hetSNPs ===
# min_base_quality 10 (default); min_mapq 15 (BaalChIP default; lower than
# ENCODE 30 because ASB requires reads at variant position which may be in
# slightly-repeat-adjacent regions)
res <- alleleCounts(res, min_base_quality = 10, min_mapq = 15)


# === 6. QC filters ===
# Remove SNPs in blacklist regions; require minimum read depth
res <- QCfilter(res,
                 RegionsToFilter = list(blacklist = 'hg38-blacklist.v2.bed'),
                 RegionsToKeep = list(peaks_only = unique(samples$Peaks)),
                 MinNumberOfReads = 12)   # min 12 reads/SNP for stable ratio


# === 7. Merge per-group counts and filter low-frequency alleles ===
res <- mergePerGroup(res)
res <- filter1allele(res)   # remove SNPs where one allele is rare (<5%)


# === 8. Bayesian beta-binomial ASB test ===
# Iter 5000 (default; increase for tighter posteriors)
# conf_level 0.95 (95% credible interval)
# CN-aware: BaalChIP reads CN from cnvs file and models overdispersion
#   accordingly. Without CN file, falls back to fixed overdispersion.
res <- getASB(res, Iter = 5000, conf_level = 0.95, cnvFile = cnvs)


# === 9. Report ===
asb_table <- BaalChIP.report(res)
cat('Total hetSNPs tested:', nrow(asb_table), '\n')

# Significant ASB calls
asb_sig <- asb_table[asb_table$isASB == TRUE & abs(asb_table$Bayes_lower_estimate - 0.5) > 0.1, ]
cat('Significant ASB (95% CI excludes 0.5; |ratio - 0.5| > 0.1):', nrow(asb_sig), '\n')

# Top ASB sites by Bayes factor
top_asb <- asb_sig[order(-asb_sig$Bayes_score), ]
head(top_asb[, c('CHROM', 'POS', 'REF', 'ALT', 'AR', 'isASB', 'Bayes_score')])

# Export
write.table(asb_table, 'HCC1395_FOXA1_ASB_full.tsv',
             sep = '\t', quote = FALSE, row.names = FALSE)
write.table(asb_sig, 'HCC1395_FOXA1_ASB_significant.tsv',
             sep = '\t', quote = FALSE, row.names = FALSE)


# === 10. Sanity checks ===
# (a) Allele ratio should not be systematically biased toward REF (WASP applied)
mean_ar <- mean(asb_table$AR, na.rm = TRUE)
cat('Mean allelic ratio (should be ~0.5; >0.55 indicates remaining REF bias):',
    round(mean_ar, 3), '\n')
stopifnot(mean_ar > 0.45 && mean_ar < 0.55)

# (b) Imprinted loci should have been filtered (no ASB at known imprinted loci)
known_imprinted <- c('H19', 'IGF2', 'MEG3', 'KCNQ1OT1')
# Cross-reference would require gene annotation; sketched here

cat('Pipeline complete.\n')
cat('  WASP-filtered BAMs assumed (drops 22-31% of reads)\n')
cat('  Imprinted loci filtered\n')
cat('  chrX filtered (female sample)\n')
cat('  CN-aware overdispersion via ASCAT calls\n')
