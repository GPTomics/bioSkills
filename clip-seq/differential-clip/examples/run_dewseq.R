# Reference: DEWSeq 1.18+, DESeq2 1.44+ | Verify API if version differs
# Differential CLIP-seq with DEWSeq (window-level NB GLM with SMInput interaction).
# Pre-requisite: htseq-clip count matrices generated upstream from dedup BAMs.

library(DEWSeq)
library(DESeq2)

# htseq-clip outputs - already merged into a count matrix
counts <- read.table('merged_counts.tsv', sep='\t', header=TRUE, row.names=1, check.names=FALSE)
annot <- 'annotation_windows.bed'   # from htseq-clip extract

# Sample metadata - each sample is one BAM
# type: ip vs sminput; condition: treatment vs control
colData <- data.frame(
    row.names = colnames(counts),
    type = c('ip','ip','ip','ip','sminput','sminput','sminput','sminput'),
    condition = c('treat','treat','ctrl','ctrl','treat','treat','ctrl','ctrl'),
    replicate = c('r1','r2','r1','r2','r1','r2','r1','r2')
)

# Build DEWSeq dataset
# The interaction design `~ type + condition + type:condition` tests whether
# the IP/SMInput ratio shifts across conditions - this is the biologically
# meaningful test for differential binding.
# Naive `~ condition` would test count differences in any sample, confounding
# binding changes with expression changes.
dds <- DESeqDataSetFromSlidingWindows(
    countData = counts,
    colData = colData,
    annotObj = annot,
    design = ~ type + condition + type:condition
)

# Pre-filter low-count windows (DEWSeq vignette recommendation)
keep <- rowSums(counts(dds) >= 5) >= 4    # >= 5 reads in >= 4 samples
dds <- dds[keep, ]

# Run DESeq2 within DEWSeq
dds <- DESeq(dds)

# Extract the interaction term - this is "does IP-vs-input differ in treat vs ctrl?"
res <- results(dds, name = 'typeip.conditiontreat')

# Standard CLIP differential thresholds: padj < 0.05, |log2FC| > 1
sig <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
cat('Significant windows:', nrow(sig), '\n')

# Convert window IDs back to BED coords
sig_df <- as.data.frame(sig)
sig_df$window_id <- rownames(sig_df)
# DEWSeq window naming: chr_start_end (verify your htseq-clip output format)
parts <- strsplit(sig_df$window_id, '_')
sig_df$chr <- sapply(parts, function(x) x[1])
sig_df$start <- as.numeric(sapply(parts, function(x) x[length(x)-1]))
sig_df$end <- as.numeric(sapply(parts, function(x) x[length(x)]))

# Write significant-window BED
out <- sig_df[, c('chr','start','end','window_id','log2FoldChange','padj')]
write.table(out, 'sig_windows.bed', sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# IMPORTANT: aggregate adjacent significant windows into biological regions
# Do this with bedtools merge AFTER writing:
#   sort -k1,1 -k2,2n sig_windows.bed | bedtools merge -d 100 -c 5,6 -o mean,min > sig_regions.bed
cat('\nNext: aggregate adjacent significant windows into regions:\n')
cat('  sort -k1,1 -k2,2n sig_windows.bed | bedtools merge -d 100 -c 5,6 -o mean,min > sig_regions.bed\n')

# Sanity check: MA plot for diagnostics
png('dewseq_MA.png', width=800, height=600)
plotMA(res, alpha=0.05, ylim=c(-5, 5), main='DEWSeq differential binding')
dev.off()

# Sanity check: top differential windows
top <- sig[order(sig$padj), ][1:20, ]
print(top)

cat('\nDone. See sig_windows.bed and dewseq_MA.png\n')
cat('Validate with: motif analysis on differential regions; orthogonal method (edgeR peak-level on same data)\n')
