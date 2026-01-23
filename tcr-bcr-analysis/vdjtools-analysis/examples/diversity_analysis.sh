#!/bin/bash
# VDJtools diversity and comparison analysis

METADATA=$1
OUTPUT_DIR=${2:-"vdjtools_output"}
VDJTOOLS_JAR=${VDJTOOLS_JAR:-"vdjtools-1.2.1/vdjtools-1.2.1.jar"}

mkdir -p $OUTPUT_DIR

echo "VDJtools analysis"
echo "Metadata: $METADATA"

# Calculate diversity statistics
# Outputs Shannon, Simpson, Chao1, Gini, d50, etc.
echo "Calculating diversity metrics..."
java -Xmx4g -jar $VDJTOOLS_JAR CalcDiversityStats \
    -m $METADATA \
    $OUTPUT_DIR/diversity

# Diversity metrics interpretation:
# - Shannon: Higher = more diverse (typical range 5-15 for healthy repertoire)
# - Gini: 0 = perfectly equal, 1 = one clone dominates (healthy: 0.3-0.6)
# - d50: Number of clones comprising 50% of repertoire (lower = oligoclonal)
# - Chao1: Estimated total richness including unseen clones

# Calculate spectratype (CDR3 length distribution)
# Normal repertoire shows Gaussian distribution centered at ~45 nt
echo "Calculating spectratype..."
java -Xmx4g -jar $VDJTOOLS_JAR CalcSpectratype \
    -m $METADATA \
    $OUTPUT_DIR/spectratype

# Calculate V and J gene usage
echo "Calculating segment usage..."
java -Xmx4g -jar $VDJTOOLS_JAR CalcSegmentUsage \
    -m $METADATA \
    $OUTPUT_DIR/segments

# Calculate pairwise overlap between samples
# Uses amino acid CDR3 for matching
echo "Calculating pairwise distances..."
java -Xmx4g -jar $VDJTOOLS_JAR CalcPairwiseDistances \
    -m $METADATA \
    -i aa \
    $OUTPUT_DIR/overlap

# Overlap metrics:
# - F2: Frequency-weighted Jaccard (recommended)
# - Jaccard: Fraction of shared unique clonotypes
# - MorisitaHorn: Accounts for clone abundances

# Find public clones (shared across multiple samples)
echo "Finding shared clonotypes..."
java -Xmx4g -jar $VDJTOOLS_JAR JoinSamples \
    -m $METADATA \
    -p \
    $OUTPUT_DIR/shared

echo ""
echo "Analysis complete!"
echo "Results in $OUTPUT_DIR:"
ls -lh $OUTPUT_DIR/

echo ""
echo "Key output files:"
echo "  diversity.txt - Diversity metrics per sample"
echo "  spectratype.txt - CDR3 length distributions"
echo "  segments.txt - V/J gene usage"
echo "  overlap.intersect.batch.txt - Pairwise overlap matrix"
echo "  shared.*.txt - Shared clonotypes"
