library(methylKit)
library(annotatr)
library(GenomicRanges)

file_list <- list('ctrl1.bismark.cov.gz', 'ctrl2.bismark.cov.gz',
                   'treat1.bismark.cov.gz', 'treat2.bismark.cov.gz')
sample_ids <- c('ctrl_1', 'ctrl_2', 'treat_1', 'treat_2')
treatment <- c(0, 0, 1, 1)

meth_obj <- methRead(location = file_list, sample.id = as.list(sample_ids), treatment = treatment,
                      assembly = 'hg38', pipeline = 'bismarkCoverage')

meth_filt <- filterByCoverage(meth_obj, lo.count = 10, hi.perc = 99.9)

tiles <- tileMethylCounts(meth_filt, win.size = 1000, step.size = 1000, cov.bases = 3)

tiles_united <- unite(tiles, destrand = TRUE)

diff_tiles <- calculateDiffMeth(tiles_united, overdispersion = 'MN', mc.cores = 4)

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
