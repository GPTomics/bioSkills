# Reference: inferCNV 1.18+ | Verify API if version differs
# Reference-based CNV inference to separate malignant from normal cells.
# Expects: a raw genes-by-cells counts matrix, a cell-to-group annotation file
# (cell <tab> group), and a gene-ordering file (gene <tab> chr <tab> start <tab> stop).
# The normal/reference groups must be confident in-sample non-malignant lineages
# (e.g. T cells, myeloid) identified by markers, never tumor-contaminated cells.
library(infercnv)

# cutoff = 0.1 suits 10x/droplet (sparse); use 1 for full-length Smart-seq.
# ref_group_names names the copy-neutral baseline groups from the annotation file.
infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = 'counts.matrix',
    annotations_file = 'cell_annotations.txt',
    delim = '\t',
    gene_order_file = 'gene_ordering.txt',
    ref_group_names = c('Tcell', 'Myeloid'))

# denoise improves signal-to-noise; HMM predicts discrete copy states per region.
# HMM_type 'i6' (default) = six copy states; 'i3' = del/neutral/amp, more robust.
# out_dir collects the heatmap and per-region state calls.
infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = 'infercnv_out',
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = TRUE,
    HMM_type = 'i6',
    num_threads = 4)

# inferCNV gives no automatic per-cell class, unlike copyKAT/SCEVAN. Derive one by
# scoring each observation cell by squared deviation of its denoised profile from
# copy-neutral, then threshold (or rerun with analysis_mode = 'subclusters' to split
# observation cells on CNV signal). add_to_seurat writes the CNA metadata onto a
# Seurat object for plotting.
obs <- read.table('infercnv_out/infercnv.observations.txt', header = TRUE, row.names = 1)
cnv_score <- colSums((obs - 1)^2)
malignant <- cnv_score > quantile(cnv_score, 0.5)

# Malignant calling is a clustering decision on a noisy proxy: validate the
# observation-vs-reference split against lineage markers or mutations, and
# treat any subclones in the heatmap as hypotheses, not measurements.
message('inferCNV complete; inspect infercnv_out/ heatmap, HMM states, and per-cell cnv_score')
