# Reference: DESeq2 1.40+, tximport 1.28+ | Verify API if version differs
# Invoked by the Snakefile `script:` directive. Snakemake injects an S4 `snakemake`
# object exposing input/output/log by name - no argument parsing needed.

library(tximport)
library(DESeq2)

log_con <- file(snakemake@log[[1]], open = 'wt')
sink(log_con)
sink(log_con, type = 'message')

quant_dirs <- snakemake@input[['quants']]
files <- file.path(quant_dirs, 'quant.sf')
names(files) <- basename(quant_dirs)

coldata <- read.csv(snakemake@input[['metadata']], row.names = 1)
coldata <- coldata[names(files), , drop = FALSE]   # align metadata rows to the salmon dirs

# txOut=TRUE keeps transcript-level estimates; pass tx2gene to summarize to genes in a real run.
txi <- tximport(files, type = 'salmon', txOut = TRUE)

# ~condition assumes a 'condition' column in the metadata; set the reference level as needed.
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~condition)
dds <- DESeq(dds)
res <- results(dds)

write.csv(as.data.frame(res[order(res$padj), ]), snakemake@output[['results']])
write.csv(counts(dds, normalized = TRUE), snakemake@output[['normalized']])
