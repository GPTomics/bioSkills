#!/bin/bash
# Reference: VDJtools 1.2.x | Verify API if version differs
# Depth-normalized VDJtools workflow: convert -> filter -> downsample -> diversity/overlap.
# The DownSample step is not optional: diversity, clonality and overlap are all depth-dependent,
# so any cross-sample comparison on un-normalized libraries measures depth, not biology.

set -euo pipefail

METADATA=${1:?"usage: diversity_analysis.sh metadata.txt [output_dir] [downsample_reads]"}
OUTPUT_DIR=${2:-vdjtools_output}
# Downsample target = the common depth all samples are reduced to. Set at/below the smallest
# library's read count; 1e5 is a typical bulk-TCR floor. Read the actual per-sample depths first.
DOWNSAMPLE_READS=${3:-100000}
VDJTOOLS_JAR=${VDJTOOLS_JAR:-vdjtools-1.2.1/vdjtools-1.2.1.jar}
JAVA_MEM=${JAVA_MEM:-4g}

mkdir -p "$OUTPUT_DIR"
run() { java -Xmx"$JAVA_MEM" -jar "$VDJTOOLS_JAR" "$@"; }

echo "[1/6] Keeping functional clonotypes only (in-frame, no stop codon)"
run FilterNonFunctional -m "$METADATA" "$OUTPUT_DIR/functional/"

echo "[2/6] Downsampling every sample to $DOWNSAMPLE_READS reads (common depth)"
# Capital S; -x/--size is the target read count. This is what makes samples comparable.
run DownSample -x "$DOWNSAMPLE_READS" -m "$METADATA" "$OUTPUT_DIR/downsampled/"

# Use the downsampled metadata for all downstream comparisons.
DS_META="$OUTPUT_DIR/downsampled/metadata.txt"

echo "[3/6] Diversity statistics (report the RESAMPLED table for cross-sample claims)"
# Emits diversity.<i>.txt (original) and diversity.<i>.resampled.txt (depth-normalized).
# Columns: observedDiversity, chao1, chaoE, efronThisted, shannonWienerIndex,
# normalizedShannonWienerIndex, inverseSimpson, d50 (each _mean/_std). Report a Hill profile:
# q=0 chaoE/observedDiversity, q=1 shannonWienerIndex, q=2 inverseSimpson. Clonality = 1 - normalizedShannonWienerIndex.
run CalcDiversityStats -m "$DS_META" "$OUTPUT_DIR/diversity"

echo "[4/6] Rarefaction curves (compare diversity at a common x, never at curve endpoints)"
run RarefactionPlot -m "$DS_META" "$OUTPUT_DIR/rarefaction"

echo "[5/6] CDR3-length spectratype and V/J segment usage"
# Spectratype: Gaussian = polyclonal; skew/spikes = clonal expansion. Compare V/J usage only within one protocol.
run CalcSpectratype -m "$DS_META" "$OUTPUT_DIR/spectratype"
run CalcSegmentUsage -m "$DS_META" "$OUTPUT_DIR/segments"

echo "[6/6] Pairwise overlap (fix the -i match key; report depth-robust MorisitaHorn/F2, not Jaccard)"
# -i aa keys on CDR3 amino acid; use nt/ntVJ for stricter, convergence-resistant sharing. Hold constant study-wide.
run CalcPairwiseDistances -i aa -m "$DS_META" "$OUTPUT_DIR/overlap"
run ClusterSamples -e MorisitaHorn "$OUTPUT_DIR/overlap" "$OUTPUT_DIR/clustered"

echo "Done. Diversity: use $OUTPUT_DIR/diversity.*.resampled.txt for comparisons."
ls -lh "$OUTPUT_DIR"
