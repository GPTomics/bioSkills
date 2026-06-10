# Reference: methylKit 1.28+, annotatr 1.28+ | Verify API if version differs
# methylKit fixed-tile DMR SCREEN. Fast and reproducible, but the per-tile q is a
# SCREENING q (fixed windows sidestep boundary selection, yet inter-tile correlation is
# ignored) - NOT a selection-corrected region FDR. For the headline call use dmrseq
# (see dmr_dmrseq.R). Thresholds below are CONVENTIONS to report, not derived truths.
library(methylKit)
library(annotatr)
library(GenomicRanges)
library(rtracklayer)

# Public methylation data sources:
# - GEO GSE86833 (WGBS), GSE105018 (RRBS tumor vs normal)
# - Bioconductor bsseqData package ships example .cov files

file_list <- list('ctrl1.bismark.cov.gz', 'ctrl2.bismark.cov.gz',
                   'treat1.bismark.cov.gz', 'treat2.bismark.cov.gz')
sample_ids <- as.list(c('ctrl_1', 'ctrl_2', 'treat_1', 'treat_2'))
treatment <- c(0, 0, 1, 1)

meth_obj <- methRead(location = file_list, sample.id = sample_ids, treatment = treatment,
                     assembly = 'hg38', pipeline = 'bismarkCoverage')

# Filter BEFORE tiling so low-coverage and library-size artifacts do not propagate.
lo_count <- 10     # single-CpG beta is a coin flip below ~10x coverage
hi_perc <- 99.9    # drop the top 0.1% coverage (PCR/repeat artifacts)
meth_filt <- filterByCoverage(meth_obj, lo.count = lo_count, hi.perc = hi_perc)
meth_norm <- normalizeCoverage(meth_filt, method = 'median')

win_size <- 1000   # 1kb windows; use 500bp for RRBS / finer resolution
cov_bases <- 3     # require >=3 CpGs per tile; the DEFAULT 0 lets single-CpG tiles through
tiles <- tileMethylCounts(meth_norm, win.size = win_size, step.size = win_size, cov.bases = cov_bases)

tiles_united <- unite(tiles, destrand = FALSE)   # destrand is for single CpGs with strand info, not tiled regions (no-op here)

# overdispersion='MN' is recommended for replicated data; note MN SILENTLY forces the
# F-test, so test='Chisq' would be ignored. adjust='BH' for cross-tool comparability
# (methylKit's default adjust='SLIM' is its own q-method, not Benjamini-Hochberg).
diff_tiles <- calculateDiffMeth(tiles_united, overdispersion = 'MN', adjust = 'BH', mc.cores = 4)

difference <- 25   # delta-beta floor (%); a CONVENTION, not derived - report and justify
qvalue <- 0.01     # FDR floor across the tile multiple-testing burden
dmrs <- getMethylDiff(diff_tiles, difference = difference, qvalue = qvalue)
dmrs_hyper <- getMethylDiff(diff_tiles, difference = difference, qvalue = qvalue, type = 'hyper')
dmrs_hypo <- getMethylDiff(diff_tiles, difference = difference, qvalue = qvalue, type = 'hypo')

sprintf('Tiles tested: %d | DMRs: %d (hyper %d, hypo %d) - SCREENING q, corroborate with dmrseq',
        nrow(diff_tiles), nrow(dmrs), nrow(dmrs_hyper), nrow(dmrs_hypo))

dmr_gr <- as(dmrs, 'GRanges')
annots <- build_annotations(genome = 'hg38', annotations = c('hg38_basicgenes', 'hg38_cpg_islands'))
dmr_annotated <- annotate_regions(regions = dmr_gr, annotations = annots, ignore.strand = TRUE)
# annotate_regions returns one row per (DMR, feature) overlap; collapse for gene-level summaries.
# For GO enrichment use missMethyl::goregion (CpG-density-bias aware), NOT a plain hypergeometric test.

dmr_annotated
