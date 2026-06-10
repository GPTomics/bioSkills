#!/bin/bash
# Reference: QIIME2 2026.1+ | Verify release/API if version differs
#
# Artifact-lifecycle walk-through: import -> orchestrate -> view/replay -> export.
# This script demonstrates the FRAMEWORK mechanics QIIME2 owns. Every scientific
# parameter choice is DEFERRED to the owning sibling skill (noted inline) and the
# values here are placeholders, not recommendations.
#
# Release model is a moving target: QIIME2 is calendar-versioned, ships as separate
# distributions (amplicon/moshpit/pathogenome/tiny), the framework is now `rachis`
# (2026.1), and `amplicon` is renamed `qiime2` in 2026.4. Run `qiime info` and verify
# distribution/action names against the current docs before pinning.

set -euo pipefail

MANIFEST='manifest.tsv'        # V2 TSV: sample-id<TAB>forward-absolute-filepath<TAB>reverse-absolute-filepath
METADATA='metadata.tsv'        # ID column + a #q2:types row annotating integer ID/batch columns categorical
CLASSIFIER='classifier.qza'    # MUST be trained for THIS release; data.qiime2.org/<release>/common/... is release-namespaced
SAMPLING_DEPTH=10000           # placeholder; pick the real depth from alpha-rarefaction -> diversity-analysis
OUT='qiime2_results'

mkdir -p "$OUT"

# 1. IMPORT - typed, provenance-rooted artifact. Phred offset is baked into the format name
#    (Phred33V2 = modern Illumina). EMP-multiplexed data needs EMPPairedEndSequences + qiime demux instead.
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MANIFEST" \
    --input-format PairedEndFastqManifestPhred33V2 \
    --output-path "$OUT/demux.qza"

# Per-base quality summary that drives the truncation choice (the choice belongs to amplicon-processing)
qiime demux summarize --i-data "$OUT/demux.qza" --o-visualization "$OUT/demux.qzv"

# 2. DENOISE - trunc/trim/maxEE and DADA2-vs-Deblur are DEFERRED to amplicon-processing
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "$OUT/demux.qza" \
    --p-trunc-len-f 0 --p-trunc-len-r 0 \
    --o-table "$OUT/table.qza" \
    --o-representative-sequences "$OUT/rep-seqs.qza" \
    --o-denoising-stats "$OUT/stats.qza"

# 3. PEEK / VALIDATE - confirm the semantic types before wiring downstream actions.
#    A type error downstream is the framework guard working; fix the upstream action, do not re-import.
qiime tools peek "$OUT/table.qza"                          # expect Type: FeatureTable[Frequency]
qiime tools validate "$OUT/table.qza" --level max          # archive integrity + format conformance

# 4. TAXONOMY - classifier and reference database are DEFERRED to taxonomy-assignment.
#    A classifier trained under another release breaks (it is pinned to its scikit-learn version).
qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads "$OUT/rep-seqs.qza" \
    --o-classification "$OUT/taxonomy.qza"

# 5. PHYLOGENY (Pipeline) - rooted tree for UniFrac / Faith PD
qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "$OUT/rep-seqs.qza" \
    --o-alignment "$OUT/aln.qza" \
    --o-masked-alignment "$OUT/masked-aln.qza" \
    --o-tree "$OUT/unrooted-tree.qza" \
    --o-rooted-tree "$OUT/rooted-tree.qza"

# 6. DIVERSITY (Pipeline) - sampling depth, metric, and rarefy-or-not are DEFERRED to diversity-analysis.
#    The PERMANOVA location-vs-dispersion (betadisper) confound also lives in diversity-analysis.
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "$OUT/rooted-tree.qza" \
    --i-table "$OUT/table.qza" \
    --p-sampling-depth "$SAMPLING_DEPTH" \
    --m-metadata-file "$METADATA" \
    --output-dir "$OUT/core-metrics"

# 7. DIFFERENTIAL ABUNDANCE - MODERN q2-composition ANCOM-BC (NOT the legacy add-pseudocount + ancom idiom).
#    Tool choice and the run-several-take-consensus message are DEFERRED to differential-abundance.
qiime composition ancombc \
    --i-table "$OUT/table.qza" \
    --m-metadata-file "$METADATA" \
    --p-formula 'group' \
    --o-differentials "$OUT/ancombc.qza"
qiime composition da-barplot --i-data "$OUT/ancombc.qza" --o-visualization "$OUT/ancombc-barplot.qzv"

# 8. PROVENANCE REPLAY - regenerate the executable commands + citations from an artifact alone.
#    Verify flag spelling with `qiime tools replay-provenance --help` (the interface is still maturing).
qiime tools replay-provenance --in-fp "$OUT/core-metrics" --out-fp "$OUT/replay.sh" --usage-driver cli
qiime tools replay-citations  --in-fp "$OUT/core-metrics" --out-fp "$OUT/citations.bib"

# 9. EXPORT - the one-way door. This DROPS the QIIME2 wrapper and the provenance; the exported TSV has no
#    history back to the reads. Export at the LAST step. For R, qiime2R::qza_to_phyloseq() reads .qza directly.
qiime tools export --input-path "$OUT/table.qza" --output-path "$OUT/exported"
biom convert -i "$OUT/exported/feature-table.biom" -o "$OUT/exported/feature-table.tsv" --to-tsv

echo "Done. View .qzv files and provenance at https://view.qiime2.org/ (no install needed)."
