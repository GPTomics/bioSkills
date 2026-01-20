#!/bin/bash
# Complete QIIME2 16S analysis workflow

set -e

# Configuration
MANIFEST="manifest.tsv"
METADATA="metadata.tsv"
CLASSIFIER="silva-138-99-nb-classifier.qza"
SAMPLING_DEPTH=10000
OUTPUT_DIR="qiime2_results"

mkdir -p "$OUTPUT_DIR"

# 1. Import data
echo "Importing data..."
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MANIFEST" \
    --output-path "$OUTPUT_DIR/demux.qza" \
    --input-format PairedEndFastqManifestPhred33V2

# 2. Denoise with DADA2
echo "Denoising with DADA2..."
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "$OUTPUT_DIR/demux.qza" \
    --p-trunc-len-f 240 \
    --p-trunc-len-r 160 \
    --p-n-threads 0 \
    --o-table "$OUTPUT_DIR/table.qza" \
    --o-representative-sequences "$OUTPUT_DIR/rep-seqs.qza" \
    --o-denoising-stats "$OUTPUT_DIR/denoising-stats.qza"

# 3. Assign taxonomy
echo "Assigning taxonomy..."
qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads "$OUTPUT_DIR/rep-seqs.qza" \
    --o-classification "$OUTPUT_DIR/taxonomy.qza"

# 4. Build phylogeny
echo "Building phylogenetic tree..."
qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "$OUTPUT_DIR/rep-seqs.qza" \
    --o-alignment "$OUTPUT_DIR/aligned-rep-seqs.qza" \
    --o-masked-alignment "$OUTPUT_DIR/masked-aligned-rep-seqs.qza" \
    --o-tree "$OUTPUT_DIR/unrooted-tree.qza" \
    --o-rooted-tree "$OUTPUT_DIR/rooted-tree.qza"

# 5. Diversity analysis
echo "Running diversity analysis..."
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "$OUTPUT_DIR/rooted-tree.qza" \
    --i-table "$OUTPUT_DIR/table.qza" \
    --p-sampling-depth "$SAMPLING_DEPTH" \
    --m-metadata-file "$METADATA" \
    --output-dir "$OUTPUT_DIR/diversity"

# 6. Alpha diversity significance
qiime diversity alpha-group-significance \
    --i-alpha-diversity "$OUTPUT_DIR/diversity/shannon_vector.qza" \
    --m-metadata-file "$METADATA" \
    --o-visualization "$OUTPUT_DIR/shannon-significance.qzv"

# 7. Beta diversity significance
qiime diversity beta-group-significance \
    --i-distance-matrix "$OUTPUT_DIR/diversity/weighted_unifrac_distance_matrix.qza" \
    --m-metadata-file "$METADATA" \
    --m-metadata-column Group \
    --p-method permanova \
    --o-visualization "$OUTPUT_DIR/permanova.qzv"

# 8. Taxonomic barplot
qiime taxa barplot \
    --i-table "$OUTPUT_DIR/table.qza" \
    --i-taxonomy "$OUTPUT_DIR/taxonomy.qza" \
    --m-metadata-file "$METADATA" \
    --o-visualization "$OUTPUT_DIR/taxa-barplot.qzv"

echo "Analysis complete! View results at https://view.qiime2.org/"
