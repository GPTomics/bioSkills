#!/usr/bin/env bash
# Reference: QIIME2 2024.2+ | Verify API if version differs
# Amplicon diversity the QIIME2 way: pick a sampling depth from the evidence, build a defensible
# tree, then run core-metrics-phylogenetic. The sampling depth SILENTLY DROPS samples below it.
set -euo pipefail

TABLE=table.qza            # FeatureTable[Frequency] from DADA2/Deblur
REPSEQS=rep-seqs.qza       # FeatureData[Sequence]
METADATA=metadata.tsv
SEPP_REF=sepp-refs-gg-13-8.qza   # full-length reference package for fragment insertion
THREADS=4

# --- Knob 1: read the per-sample frequencies and the rarefaction plateau to CHOOSE the depth ---
qiime feature-table summarize --i-table "$TABLE" --o-visualization table.qzv
qiime diversity alpha-rarefaction \
    --i-table "$TABLE" \
    --p-max-depth 20000 \
    --m-metadata-file "$METADATA" \
    --o-visualization alpha-rarefaction.qzv
# Inspect table.qzv (per-sample counts) and alpha-rarefaction.qzv (the curve + survivors-per-depth panel),
# then set SAMPLING_DEPTH on the observed-features plateau that retains an acceptable fraction of samples.
SAMPLING_DEPTH=10000

# --- Knob 2: the tree. SEPP fragment-insertion into a full-length reference beats a de novo build. ---
qiime fragment-insertion sepp \
    --i-representative-sequences "$REPSEQS" \
    --i-reference-database "$SEPP_REF" \
    --p-threads "$THREADS" \
    --o-tree insertion-tree.qza \
    --o-placements insertion-placements.qza
# Fragments that failed to insert MUST be removed before diversity (a second silent table-shrink):
qiime fragment-insertion filter-features \
    --i-table "$TABLE" \
    --i-tree insertion-tree.qza \
    --o-filtered-table table-sepp.qza \
    --o-removed-table removed-table.qza

# --- core-metrics: rarefies to SAMPLING_DEPTH (dropping samples below it), computes 4 alpha + 4 beta ---
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny insertion-tree.qza \
    --i-table table-sepp.qza \
    --p-sampling-depth "$SAMPLING_DEPTH" \
    --m-metadata-file "$METADATA" \
    --output-dir core-metrics-results

# --- Knob 3: test BOTH UniFrac variants for location, and run permdisp for dispersion ---
for METRIC in weighted_unifrac unweighted_unifrac; do
    qiime diversity beta-group-significance \
        --i-distance-matrix "core-metrics-results/${METRIC}_distance_matrix.qza" \
        --m-metadata-file "$METADATA" \
        --m-metadata-column Group \
        --p-method permanova --p-pairwise \
        --o-visualization "${METRIC}-permanova.qzv"
    # MANDATORY dispersion check alongside PERMANOVA:
    qiime diversity beta-group-significance \
        --i-distance-matrix "core-metrics-results/${METRIC}_distance_matrix.qza" \
        --m-metadata-file "$METADATA" \
        --m-metadata-column Group \
        --p-method permdisp \
        --o-visualization "${METRIC}-permdisp.qzv"
done

echo 'Done. Report the sampling depth, the dropped-sample count, and both PERMANOVA + permdisp results.'
