# Reference: methylKit 1.28+, GenomicRanges 1.54+ | Verify API if version differs
# The methylKit import-to-results spine: import -> filter -> normalize -> unite ->
# QC -> calculateDiffMeth -> getMethylDiff, for both per-CpG (DMC) and tile (DMR) results.
library(methylKit)

file_list <- list('ctrl1.cov.gz', 'ctrl2.cov.gz', 'treat1.cov.gz', 'treat2.cov.gz')
sample_ids <- list('ctrl_1', 'ctrl_2', 'treat_1', 'treat_2')
treatment <- c(0, 0, 1, 1)

# pipeline='bismarkCoverage' reads .cov (no strand); use 'bismarkCytosineReport'
# when destranding is needed (the report carries strand + context).
meth_obj <- methRead(file_list, sample.id=sample_ids, treatment=treatment,
                     assembly='hg38', context='CpG', pipeline='bismarkCoverage')

getMethylationStats(meth_obj[[1]], plot=FALSE, both.strands=FALSE)
getCoverageStats(meth_obj[[1]], plot=FALSE, both.strands=FALSE)

lo_count <- 10      # below ~10x a single-CpG percentage is a coin flip
hi_perc <- 99.9     # drop the top 0.1% coverage (PCR/repeat artifacts)
meth_filt <- filterByCoverage(meth_obj, lo.count=lo_count, lo.perc=NULL, hi.count=NULL, hi.perc=hi_perc)

# normalize BEFORE unite/tile so a deeper sample does not look more confident
meth_norm <- normalizeCoverage(meth_filt, method='median')

# destrand=FALSE here because .cov carries no strand; use TRUE only with a cytosine report
meth_united <- unite(meth_norm, destrand=FALSE)

getCorrelation(meth_united, plot=FALSE)
clusterSamples(meth_united, dist='correlation', method='ward.D', plot=FALSE)

# overdispersion='MN' adds the beta-binomial layer and FORCES the F-test (test= is ignored);
# adjust='BH' for cross-tool comparability (default is SLIM, methylKit-specific).
diff_meth <- calculateDiffMeth(meth_united, overdispersion='MN', adjust='BH', mc.cores=1)

difference <- 25    # 25% is tutorial convention, NOT derived; justify and report per study
qvalue <- 0.01      # FDR floor across the per-CpG multiple-testing burden
dmcs <- getMethylDiff(diff_meth, difference=difference, qvalue=qvalue)
dmcs_hyper <- getMethylDiff(diff_meth, difference=difference, qvalue=qvalue, type='hyper')
dmcs_hypo <- getMethylDiff(diff_meth, difference=difference, qvalue=qvalue, type='hypo')

nrow(dmcs)
nrow(dmcs_hyper)
nrow(dmcs_hypo)

# Tile screen: tile the FILTERED/NORMALIZED object, require real CpG support (cov.bases>=3).
# The per-tile q is a screen, NOT a selection-corrected region FDR (see dmr-detection).
tiles <- tileMethylCounts(meth_norm, win.size=1000, step.size=1000, cov.bases=3)
tiles_united <- unite(tiles, destrand=FALSE)
diff_tiles <- calculateDiffMeth(tiles_united, overdispersion='MN', adjust='BH', mc.cores=1)
dmrs <- getMethylDiff(diff_tiles, difference=difference, qvalue=qvalue)
nrow(dmrs)

# pool() sums replicates into one pseudo-sample per group; it DESTROYS biological
# replication and must NOT feed the reported test - exploratory visualization only.
# meth_pooled <- pool(meth_united, sample.ids=c('control', 'treatment'))

out_file <- tempfile(fileext='.csv')
write.csv(getData(dmcs), out_file, row.names=FALSE)
