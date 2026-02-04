#!/usr/bin/env Rscript
# Differential splicing analysis with leafcutter
# Intron cluster-based approach for novel junction discovery

library(leafcutter)

# Step 1: Convert BAMs to junction files using regtools
# Run from command line:
# for bam in *.bam; do
#     regtools junctions extract -a 8 -m 50 -s 0 $bam -o ${bam%.bam}.junc
# done

# Step 2: Create junction file list
junc_files <- list.files(pattern = '\\.junc$', full.names = TRUE)
writeLines(junc_files, 'juncfiles.txt')

# Step 3: Cluster introns (run from command line)
# python leafcutter_cluster_regtools.py \
#     -j juncfiles.txt \
#     -m 50 \       # Min reads per cluster
#     -o leafcutter \
#     -l 500000     # Max intron length

# Step 4: Create groups file for differential analysis
groups <- data.frame(
    sample = c('sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6'),
    group = c('control', 'control', 'control', 'treatment', 'treatment', 'treatment')
)
write.table(groups, 'groups.txt', sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

# Step 5: Run differential splicing (command line)
# Rscript leafcutter_ds.R \
#     leafcutter_perind_numers.counts.gz \
#     groups.txt \
#     -o differential \
#     -e annotation_exons.txt.gz

# Step 6: Load and analyze results
load_leafcutter_results <- function(cluster_sig_file, effect_sizes_file) {
    sig <- read.table(cluster_sig_file, header = TRUE, sep = '\t')
    effects <- read.table(effect_sizes_file, header = TRUE, sep = '\t')

    # Filter significant clusters
    # FDR < 0.05 (p.adjust column)
    significant <- sig[sig$p.adjust < 0.05, ]

    cat(sprintf('Total clusters tested: %d\n', nrow(sig)))
    cat(sprintf('Significant clusters (FDR < 0.05): %d\n', nrow(significant)))

    # Merge with effect sizes
    results <- merge(significant, effects, by = 'cluster', all.x = TRUE)

    # Sort by significance
    results <- results[order(results$p.adjust), ]

    return(results)
}

# Example: Parse results
# results <- load_leafcutter_results(
#     'differential_cluster_significance.txt',
#     'differential_effect_sizes.txt'
# )
# head(results)

# Annotate with gene information
annotate_clusters <- function(results, exon_file) {
    exons <- read.table(exon_file, header = TRUE)

    # Match cluster coordinates to genes
    # Leafcutter cluster format: chr:start:end:clu_N
    results$chr <- sapply(strsplit(results$cluster, ':'), '[', 1)
    results$start <- as.numeric(sapply(strsplit(results$cluster, ':'), '[', 2))
    results$end <- as.numeric(sapply(strsplit(results$cluster, ':'), '[', 3))

    return(results)
}

cat('leafcutter differential splicing analysis pipeline\n')
cat('Run clustering and differential testing steps from command line\n')
