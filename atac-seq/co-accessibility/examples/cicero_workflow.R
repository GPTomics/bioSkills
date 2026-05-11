#!/usr/bin/env Rscript
# Reference: cicero 1.20+, monocle3 1.3+, BSgenome.Hsapiens.UCSC.hg38 1.4+, GenomicRanges 1.54+ | Verify API if version differs
# Cicero co-accessibility: scATAC peak matrix -> metacell aggregation -> graphical lasso ->
# peak-pair connection scores -> enhancer-gene candidate pairs.

suppressPackageStartupMessages({
    library(cicero); library(monocle3); library(GenomicRanges)
    library(BSgenome.Hsapiens.UCSC.hg38); library(rtracklayer)
})

run_cicero_pipeline <- function(peak_matrix, peak_metadata, cell_metadata,
                                output_prefix='cicero',
                                k_metacell=50, window=500000, coaccess_threshold=0.25,
                                tss_bed='gencode_v29_protein_coding_tss.bed',
                                tss_pad=2000) {

    # 1. Build CDS from inputs (peaks x cells matrix; expects integer 0/1)
    input_cds <- new_cell_data_set(peak_matrix,
                                   cell_metadata=cell_metadata,
                                   gene_metadata=peak_metadata)
    cat(sprintf('Loaded: %d peaks, %d cells\n', nrow(input_cds), ncol(input_cds)))

    # 2. Dimensionality reduction (LSI for sparse binary data; then UMAP)
    input_cds <- detect_genes(input_cds)
    input_cds <- estimate_size_factors(input_cds)
    input_cds <- preprocess_cds(input_cds, method='LSI')
    input_cds <- reduce_dimension(input_cds, reduction_method='UMAP',
                                  preprocess_method='LSI')

    # 3. Build metacells via k-NN (default k=50). Smaller k -> more variability captured but slower.
    umap_coords <- reducedDims(input_cds)$UMAP
    cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates=umap_coords, k=k_metacell)

    # 4. Run Cicero with hg38 chromosome sizes (cis only by default)
    chrs <- seqnames(seqinfo(BSgenome.Hsapiens.UCSC.hg38))
    chrs <- chrs[!grepl('_alt|_random|chrUn', chrs)]
    genome_df <- data.frame(chr=chrs, length=seqlengths(BSgenome.Hsapiens.UCSC.hg38)[chrs])

    cat('Running Cicero (this can take time on large datasets)...\n')
    conns <- run_cicero(cicero_cds, genomic_coords=genome_df,
                        window=window, sample_num=100)

    cat(sprintf('Total connections: %d\n', nrow(conns)))
    strong <- conns[conns$coaccess > coaccess_threshold & !is.na(conns$coaccess), ]
    cat(sprintf('Strong (coaccess > %.2f): %d\n', coaccess_threshold, nrow(strong)))

    # 5. Save BEDPE for visualization
    write.table(strong, sprintf('%s_connections.tsv', output_prefix),
                sep='\t', row.names=FALSE, quote=FALSE)

    # 6. Map to enhancer-gene pairs via TSS overlap
    if (file.exists(tss_bed)) {
        tss <- import(tss_bed)
        tss_extended <- resize(tss, width=2 * tss_pad, fix='center')

        # Convert connections to GRanges (each connection has two peaks)
        peak1 <- GRanges(strong$Peak1)
        peak2 <- GRanges(strong$Peak2)
        ov1 <- findOverlaps(peak1, tss_extended)
        ov2 <- findOverlaps(peak2, tss_extended)

        eg_pairs <- data.frame(
            enhancer = c(strong$Peak2[queryHits(ov1)], strong$Peak1[queryHits(ov2)]),
            gene = c(tss_extended$gene_id[subjectHits(ov1)],
                     tss_extended$gene_id[subjectHits(ov2)]),
            coaccess = c(strong$coaccess[queryHits(ov1)], strong$coaccess[queryHits(ov2)]))
        eg_pairs <- unique(eg_pairs)
        write.csv(eg_pairs, sprintf('%s_enhancer_gene_pairs.csv', output_prefix),
                  row.names=FALSE)
        cat(sprintf('Enhancer-gene pairs: %d unique\n', nrow(eg_pairs)))
    } else {
        cat(sprintf('TSS BED not found at %s; skipping enhancer-gene mapping\n', tss_bed))
    }

    invisible(list(conns=conns, strong=strong))
}

# Usage: provide peak_matrix.mtx, peak_metadata.tsv, cell_metadata.tsv as args
args <- commandArgs(trailingOnly=TRUE)
if (length(args) >= 3) {
    peak_matrix <- Matrix::readMM(args[1])
    peak_meta <- read.delim(args[2], row.names=1)
    cell_meta <- read.delim(args[3], row.names=1)
    run_cicero_pipeline(peak_matrix, peak_meta, cell_meta)
}
