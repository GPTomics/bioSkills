# Complete methylation analysis with methylKit

library(methylKit)
library(genomation)
library(ggplot2)

# Configuration
output_dir <- 'methylation_results/analysis'
dir.create(output_dir, showWarnings = FALSE)

# Sample information
files <- list(
    'methylation_results/methylation/control_1.CpG_report.txt.gz',
    'methylation_results/methylation/control_2.CpG_report.txt.gz',
    'methylation_results/methylation/treated_1.CpG_report.txt.gz',
    'methylation_results/methylation/treated_2.CpG_report.txt.gz'
)

sample_ids <- c('control_1', 'control_2', 'treated_1', 'treated_2')
treatment <- c(0, 0, 1, 1)  # 0 = control, 1 = treated

# Read Bismark cytosine reports
cat('Reading methylation data...\n')
meth_obj <- methRead(
    location = as.list(files),
    sample.id = as.list(sample_ids),
    assembly = 'hg38',
    treatment = treatment,
    context = 'CpG',
    pipeline = 'bismarkCytosineReport'
)

# QC statistics
pdf(file.path(output_dir, 'qc_plots.pdf'), width = 10, height = 8)
for (i in seq_along(meth_obj)) {
    getMethylationStats(meth_obj[[i]], plot = TRUE, both.strands = FALSE)
    getCoverageStats(meth_obj[[i]], plot = TRUE, both.strands = FALSE)
}
dev.off()

# Filter by coverage
cat('Filtering by coverage...\n')
meth_filtered <- filterByCoverage(meth_obj,
    lo.count = 10,
    lo.perc = NULL,
    hi.count = NULL,
    hi.perc = 99.9
)

# Normalize coverage
cat('Normalizing coverage...\n')
meth_norm <- normalizeCoverage(meth_filtered, method = 'median')

# Merge samples
cat('Merging samples...\n')
meth_merged <- unite(meth_norm, destrand = TRUE, min.per.group = 1L)
cat('CpGs covered in all samples:', nrow(meth_merged), '\n')

# Sample clustering
pdf(file.path(output_dir, 'sample_clustering.pdf'), width = 8, height = 6)
clusterSamples(meth_merged, dist = 'correlation', method = 'ward', plot = TRUE)
PCASamples(meth_merged)
dev.off()

# Differential methylation (per CpG)
cat('Calculating differential methylation...\n')
diff_meth <- calculateDiffMeth(meth_merged, mc.cores = 4)

# Summary
cat('\nDifferential methylation summary:\n')
print(diffMethPerChr(diff_meth, plot = FALSE))

# Significant DMCs
dmc_hyper <- getMethylDiff(diff_meth, difference = 25, qvalue = 0.01, type = 'hyper')
dmc_hypo <- getMethylDiff(diff_meth, difference = 25, qvalue = 0.01, type = 'hypo')
cat('\nHypermethylated CpGs:', nrow(dmc_hyper), '\n')
cat('Hypomethylated CpGs:', nrow(dmc_hypo), '\n')

# DMR detection using tiles
cat('\nDetecting DMRs using tiles...\n')
tiles <- tileMethylCounts(meth_merged, win.size = 1000, step.size = 1000)
diff_tiles <- calculateDiffMeth(tiles, mc.cores = 4)
dmr <- getMethylDiff(diff_tiles, difference = 25, qvalue = 0.01)
cat('Significant DMRs:', nrow(dmr), '\n')

# Export results
write.csv(as.data.frame(diff_meth), file.path(output_dir, 'all_cpg_diff.csv'), row.names = FALSE)
write.csv(as.data.frame(dmr), file.path(output_dir, 'dmr_results.csv'), row.names = FALSE)

# Export as BED
dmr_gr <- as(dmr, 'GRanges')
hyper_dmr <- dmr_gr[dmr_gr$meth.diff > 0]
hypo_dmr <- dmr_gr[dmr_gr$meth.diff < 0]

export_bed <- function(gr, filename) {
    if (length(gr) > 0) {
        df <- data.frame(
            chr = seqnames(gr),
            start = start(gr) - 1,
            end = end(gr),
            name = paste0('DMR_', seq_along(gr)),
            score = abs(gr$meth.diff)
        )
        write.table(df, filename, sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
}

export_bed(hyper_dmr, file.path(output_dir, 'hyper_dmr.bed'))
export_bed(hypo_dmr, file.path(output_dir, 'hypo_dmr.bed'))

cat('\nResults saved to:', output_dir, '\n')
