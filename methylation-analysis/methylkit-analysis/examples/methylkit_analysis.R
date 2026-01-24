library(methylKit)

file_list <- list('ctrl1.bismark.cov.gz', 'ctrl2.bismark.cov.gz',
                   'treat1.bismark.cov.gz', 'treat2.bismark.cov.gz')
sample_ids <- c('ctrl_1', 'ctrl_2', 'treat_1', 'treat_2')
treatment <- c(0, 0, 1, 1)

meth_obj <- methRead(
    location = file_list,
    sample.id = as.list(sample_ids),
    treatment = treatment,
    assembly = 'hg38',
    context = 'CpG',
    pipeline = 'bismarkCoverage'
)

getMethylationStats(meth_obj[[1]], plot = TRUE, both.strands = FALSE)
getCoverageStats(meth_obj[[1]], plot = TRUE, both.strands = FALSE)

# lo.count=10: Minimum 10 reads per CpG. Standard for reliable methylation calls.
# hi.perc=99.9: Remove top 0.1% coverage (likely PCR artifacts or repeat regions).
meth_filt <- filterByCoverage(meth_obj, lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)

meth_norm <- normalizeCoverage(meth_filt, method = 'median')

meth_united <- unite(meth_norm, destrand = TRUE)

getCorrelation(meth_united, plot = TRUE)
PCASamples(meth_united)
clusterSamples(meth_united, dist = 'correlation', method = 'ward.D', plot = TRUE)

diff_meth <- calculateDiffMeth(meth_united, overdispersion = 'MN', test = 'Chisq', mc.cores = 4)

# difference=25: Minimum 25% methylation difference. Standard threshold.
# qvalue=0.01: FDR-adjusted p-value cutoff. Use 0.05 for exploratory analysis.
dmcs <- getMethylDiff(diff_meth, difference = 25, qvalue = 0.01)
dmcs_hyper <- getMethylDiff(diff_meth, difference = 25, qvalue = 0.01, type = 'hyper')
dmcs_hypo <- getMethylDiff(diff_meth, difference = 25, qvalue = 0.01, type = 'hypo')

nrow(dmcs)
nrow(dmcs_hyper)
nrow(dmcs_hypo)

diff_df <- getData(dmcs)
write.csv(diff_df, 'dmcs_results.csv', row.names = FALSE)
