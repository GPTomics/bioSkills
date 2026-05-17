#!/usr/bin/env Rscript
# Reference: csaw 1.36+, edgeR 4.0+, rtracklayer 1.62+ | Verify API if version differs
# csaw sliding-window differential binding with composition-bias-robust
# background-bin TMM normalization. Designed for global-shift experiments
# (HDACi/BETi/EZH2i) where reads-in-peaks RLE/TMM gives wrong answers, and
# for broad marks (H3K27me3, H3K9me3) where DiffBind summit-recentering loses
# domain-level biology.

suppressPackageStartupMessages({
    library(csaw)
    library(edgeR)
    library(rtracklayer)
})

# Sample setup
bam_files <- c('ctrl_1.bam', 'ctrl_2.bam', 'ctrl_3.bam',
                'treat_1.bam', 'treat_2.bam', 'treat_3.bam')
condition <- factor(c('ctrl', 'ctrl', 'ctrl', 'treat', 'treat', 'treat'))

# Read parameters: MAPQ 30 (uniquely-mapped), paired-end, mark-and-discard
# duplicates (csaw expects unmarked dup status; pre-dedup BAM is fine), blacklist
# excluded. minq=30 matches ENCODE filter -q 30; pe='both' requires both mates.
blacklist <- import('hg38-blacklist.v2.bed')
param <- readParam(minq = 30, pe = 'both', dedup = TRUE,
                    discard = blacklist, restrict = paste0('chr', c(1:22, 'X', 'Y')))

# Window counts: narrow marks 150 bp window, 50 bp shift; broad marks 1-2 kb.
# ext=200: extend single-end reads to fragment length (use cross-correlation
# estimate); for paired-end this is ignored.
windows <- windowCounts(bam_files, width = 150, ext = 200, param = param,
                         spacing = 50)

# Filter low-abundance windows (csaw default: abundance >= 1)
filter_stat <- filterWindowsGlobal(windows, background = NULL)
keep <- filter_stat$filter > log2(3)  # 3-fold over background
windows <- windows[keep, ]

# Composition bias: TMM on 10 kb background bins (always applied for ChIP-seq)
# This is the key step that distinguishes csaw from DiffBind default. The
# background bin counts represent stable genome-wide signal; TMM size factors
# from these bins correct composition bias without assuming peaks are stable.
bg_bins <- windowCounts(bam_files, bin = TRUE, width = 10000, param = param)
windows <- normFactors(bg_bins, se.out = windows)

# OPTIONAL: trended (abundance-dependent) bias via non-linear loess.
# Use cautiously — if the trend IS the biology (e.g., treatment uniformly
# increases low-signal windows), loess removes it. Only apply when:
# - Library prep efficiency variation across abundance is suspected, OR
# - Same trend appears in IgG-only / no-IP control samples
# windows <- normOffsets(windows, type = 'loess', se.out = TRUE)

# edgeR quasi-likelihood F-test (better type-I control than Wald for small n)
y <- asDGEList(windows)
design <- model.matrix(~condition)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)
results <- glmQLFTest(fit, coef = 2)

# Merge significant adjacent windows into regions (tol=1 kb for sharp marks,
# 5 kb for broad). max.width caps regions to avoid stitching across domains.
merged <- mergeResults(windows, results$table,
                        tol = 1000, merge.args = list(max.width = 5000))

# Significant regions
sig <- merged$combined[merged$combined$FDR < 0.05, ]
cat('Significant regions (FDR < 0.05):', nrow(sig), '\n')

# Export
export(sig, 'csaw_differential_regions.bed', format = 'BED')
write.table(as.data.frame(sig), 'csaw_differential_regions.tsv',
            sep = '\t', row.names = FALSE, quote = FALSE)

# MA plot to verify normalization sanity
pdf('csaw_ma_plot.pdf')
plotMD(fit, coef = 2, status = ifelse(results$table$PValue < 0.01, 'sig', 'ns'),
        values = c('sig', 'ns'), col = c('red', 'gray'))
abline(h = 0, col = 'blue')
dev.off()
# Inspect: loess curve should be near y=0. If shifted, composition-bias
# correction failed OR there's a global shift requiring spike-in.
