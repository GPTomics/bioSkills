# Reference: GenomicFeatures 1.54+, txdbmaker 1.0+ (Bioconductor >= 3.19) | Verify API if version differs
library(GenomicFeatures)

gtf_file <- 'Homo_sapiens.GRCh38.110.gtf.gz'

cat('Creating TxDb from GTF...\n')
# makeTxDbFromGFF moved to txdbmaker in Bioconductor >= 3.19 (defunct in GenomicFeatures >= 1.61.1)
make_txdb <- if (requireNamespace('txdbmaker', quietly = TRUE)) txdbmaker::makeTxDbFromGFF else GenomicFeatures::makeTxDbFromGFF
txdb <- make_txdb(gtf_file)

cat('Extracting transcript-gene mapping...\n')
k <- keys(txdb, keytype = 'TXNAME')
tx2gene <- AnnotationDbi::select(txdb, k, c('TXNAME', 'GENEID'), 'TXNAME')
tx2gene <- tx2gene[, c('TXNAME', 'GENEID')]

cat('Saving tx2gene.csv...\n')
write.csv(tx2gene, 'tx2gene.csv', row.names = FALSE)

cat('Done! Created tx2gene.csv with', nrow(tx2gene), 'transcripts\n')
