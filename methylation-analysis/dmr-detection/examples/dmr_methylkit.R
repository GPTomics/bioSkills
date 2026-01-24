library(methylKit)
library(annotatr)
library(GenomicRanges)

# Public methylation data sources:
# - GEO: GSE86833 (WGBS from various tissues)
# - GEO: GSE105018 (RRBS tumor vs normal)
# - ENCODE: ENCSR000COQ (bismark coverage files)
# - Bioconductor: bsseqData package has example .cov files

file_list <- list('ctrl1.bismark.cov.gz', 'ctrl2.bismark.cov.gz',
                   'treat1.bismark.cov.gz', 'treat2.bismark.cov.gz')
sample_ids <- c('ctrl_1', 'ctrl_2', 'treat_1', 'treat_2')
treatment <- c(0, 0, 1, 1)

meth_obj <- methRead(location = file_list, sample.id = as.list(sample_ids), treatment = treatment,
                      assembly = 'hg38', pipeline = 'bismarkCoverage')

# lo.count=10: Minimum 10 reads per CpG. Standard for reliable methylation calls.
# hi.perc=99.9: Remove top 0.1% coverage (likely PCR duplicates or repeat regions).
meth_filt <- filterByCoverage(meth_obj, lo.count = 10, hi.perc = 99.9)

# win.size=1000: 1kb windows are standard for DMR detection. Use 500bp for finer resolution.
# cov.bases=3: Require at least 3 CpGs per window for statistical reliability.
tiles <- tileMethylCounts(meth_filt, win.size = 1000, step.size = 1000, cov.bases = 3)

tiles_united <- unite(tiles, destrand = TRUE)

diff_tiles <- calculateDiffMeth(tiles_united, overdispersion = 'MN', mc.cores = 4)

# difference=25: Minimum 25% methylation difference. Commonly used threshold.
# Use 10% for subtle effects, 30-50% for strong effects.
# qvalue=0.01: FDR-adjusted p-value cutoff. Use 0.05 for exploratory, 0.01 for stringent.
dmrs <- getMethylDiff(diff_tiles, difference = 25, qvalue = 0.01)
dmrs_hyper <- getMethylDiff(diff_tiles, difference = 25, qvalue = 0.01, type = 'hyper')
dmrs_hypo <- getMethylDiff(diff_tiles, difference = 25, qvalue = 0.01, type = 'hypo')

sprintf('Total DMRs: %d (Hyper: %d, Hypo: %d)', nrow(dmrs), nrow(dmrs_hyper), nrow(dmrs_hypo))

annots <- build_annotations(genome = 'hg38', annotations = c('hg38_basicgenes', 'hg38_cpg_islands'))
dmr_gr <- as(dmrs, 'GRanges')
dmr_annotated <- annotate_regions(regions = dmr_gr, annotations = annots, ignore.strand = TRUE)

dmr_df <- data.frame(dmr_annotated)
write.csv(dmr_df, 'dmrs_annotated.csv', row.names = FALSE)

library(rtracklayer)
export(dmr_gr, 'dmrs.bed', format = 'BED')
