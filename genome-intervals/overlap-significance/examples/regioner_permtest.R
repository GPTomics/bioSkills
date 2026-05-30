# Reference: regioneR 1.36+ (Bioconductor 3.18+) | Verify API if version differs
# Publication-grade colocalization test: permutation p + z-score with a clustering-preserving
# null, then localZScore to see where within the regions the association lives.
suppressPackageStartupMessages(library(regioneR))

N_TIMES <- 1000     # permutation count; >=1000 for a stable empirical p (Gel 2016)
LZ_WINDOW <- 10000  # bp window for localZScore shifting
LZ_STEP <- 500      # bp step for localZScore

peaks <- toGRanges('peaks.bed')
enhancers <- toGRanges('enhancers.bed')

# mask is the workspace control: exclude the ENCODE blacklist + assembly gaps from random placement
gam <- getGenomeAndMask(genome = 'hg38', mask = toGRanges('blacklist.bed'))

pt <- permTest(A = peaks, B = enhancers,
               randomize.function = circularRandomizeRegions,  # preserves clustering; use randomizeRegions if query is not autocorrelated
               evaluate.function = numOverlaps,
               genome = gam$genome, mask = gam$mask,
               ntimes = N_TIMES, count.once = TRUE)

cat(sprintf('p-value: %.4g\n', pt$numOverlaps$pval))
cat(sprintf('z-score: %.3f\n', pt$numOverlaps$zscore))

lz <- localZScore(A = peaks, B = enhancers, pt = pt, window = LZ_WINDOW, step = LZ_STEP)
lz
