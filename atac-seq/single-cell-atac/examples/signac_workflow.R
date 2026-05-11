#!/usr/bin/env Rscript
# Reference: Signac 1.13+, Seurat 5.0+, EnsDb.Hsapiens.v86 2.99+, BSgenome.Hsapiens.UCSC.hg38 1.4+ | Verify API if version differs
# Standard Signac scATAC-seq pipeline: 10X Cell Ranger output -> per-cell QC -> TF-IDF/LSI ->
# UMAP/Leiden (skipping depth-correlated dim 1) -> gene-activity scores for annotation.

suppressPackageStartupMessages({
    library(Signac); library(Seurat); library(EnsDb.Hsapiens.v86)
    library(BSgenome.Hsapiens.UCSC.hg38); library(GenomicRanges); library(ggplot2)
})

run_signac <- function(h5_file='outs/filtered_peak_bc_matrix.h5',
                       fragments_file='outs/fragments.tsv.gz',
                       metadata_file='outs/singlecell.csv',
                       output_prefix='scatac') {

    counts <- Read10X_h5(h5_file)
    metadata <- read.csv(metadata_file, header=TRUE, row.names=1)

    # Build Signac assay
    chrom_assay <- CreateChromatinAssay(
        counts=counts, sep=c(':', '-'),
        genome='hg38', fragments=fragments_file,
        annotation=GetGRangesFromEnsDb(EnsDb.Hsapiens.v86),
        min.cells=10, min.features=200)
    obj <- CreateSeuratObject(counts=chrom_assay, assay='ATAC', meta.data=metadata)

    # Per-cell QC -- looser per-cell thresholds than bulk
    obj <- NucleosomeSignal(obj)
    obj <- TSSEnrichment(obj, fast=FALSE)
    obj$pct_reads_in_peaks <- obj$peak_region_fragments / obj$passed_filters * 100
    obj$blacklist_ratio <- obj$blacklist_region_fragments / obj$peak_region_fragments

    # QC plots
    pdf(sprintf('%s_qc.pdf', output_prefix), 14, 4)
    print(VlnPlot(obj, features=c('nCount_ATAC', 'TSS.enrichment',
                                   'pct_reads_in_peaks', 'nucleosome_signal',
                                   'blacklist_ratio'), pt.size=0.1, ncol=5))
    dev.off()

    # Re-filter cellranger output at sensible per-cell thresholds (cellranger is lenient)
    obj_filt <- subset(obj,
        subset = peak_region_fragments > 1000 & peak_region_fragments < 20000 &
                 pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 &
                 nucleosome_signal < 4 & TSS.enrichment > 4)
    cat(sprintf('After QC filter: %d / %d cells\n', ncol(obj_filt), ncol(obj)))

    # Dimensionality reduction
    # CRITICAL: dims = 2:30, NOT 1:30. LSI component 1 is depth, not biology.
    obj_filt <- RunTFIDF(obj_filt)
    obj_filt <- FindTopFeatures(obj_filt, min.cutoff='q0')
    obj_filt <- RunSVD(obj_filt)

    # Confirm depth correlation: component 1 should correlate with nCount_ATAC
    depth_cor <- DepthCor(obj_filt)
    pdf(sprintf('%s_depth_cor.pdf', output_prefix), 6, 4)
    print(depth_cor)
    dev.off()
    cat('  (Component 1 should correlate ~ -1 with depth; that is why we skip it.)\n')

    obj_filt <- RunUMAP(obj_filt, reduction='lsi', dims=2:30)
    obj_filt <- FindNeighbors(obj_filt, reduction='lsi', dims=2:30)
    obj_filt <- FindClusters(obj_filt, algorithm=4, resolution=0.5)        # Leiden = algorithm 4

    pdf(sprintf('%s_umap.pdf', output_prefix), 7, 6)
    print(DimPlot(obj_filt, label=TRUE, label.size=4))
    dev.off()

    # Gene activity (approximation of expression from accessibility) for annotation
    cat('Computing gene activity scores...\n')
    gene_activities <- GeneActivity(obj_filt)
    obj_filt[['ACT']] <- CreateAssayObject(counts=gene_activities)
    DefaultAssay(obj_filt) <- 'ACT'
    obj_filt <- NormalizeData(obj_filt, normalization.method='LogNormalize',
                              scale.factor=median(obj_filt$nCount_ACT))

    saveRDS(obj_filt, sprintf('%s_signac.rds', output_prefix))
    cat(sprintf('Saved to %s_signac.rds\n', output_prefix))
    invisible(obj_filt)
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0) run_signac(args[1], args[2], args[3])
