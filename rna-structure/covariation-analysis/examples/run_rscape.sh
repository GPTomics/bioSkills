#!/bin/bash
# Reference: R-scape 2.0+ | Verify API if version differs
# Test whether an RNA alignment's secondary structure is supported by evolutionary covariation.

ALIGNMENT=$1
OUTDIR=${2:-"rscape_out"}
EVALUE=${3:-0.05}

if [ -z "$ALIGNMENT" ]; then
    echo "Usage: $0 <alignment.sto> [outdir] [evalue]"
    echo ""
    echo "Input is a Stockholm alignment. For -s the alignment must carry a #=GC SS_cons line"
    echo "(the structure to test); for --cacofold R-scape predicts a covariation-supported structure."
    exit 1
fi

mkdir -p "$OUTDIR"

echo "=== Test the given consensus structure (#=GC SS_cons) ==="
# -s scores the pairs in SS_cons against a phylogeny-aware null; E-value target default 0.05.
R-scape -s -E "$EVALUE" --outdir "$OUTDIR" "$ALIGNMENT"

echo ""
echo "=== Predict a covariation-supported structure de novo (CaCoFold) ==="
# Use when there is no trusted SS_cons to test. --cacofold is also accepted as --fold.
R-scape --cacofold -E "$EVALUE" --outdir "$OUTDIR" "$ALIGNMENT"

echo ""
echo "=== Outputs in $OUTDIR ==="
echo "  *.cov         - significantly covarying pairs (positions, score, E-value, substitutions, power)"
echo "  *.power       - per-pair statistical power analysis"
echo "  *.R2R.sto     - CaCoFold consensus structure (de novo run)"
echo "  *.svg / *.pdf - R2R diagram coloring covarying pairs on the structure"
echo ""
echo "Interpret with the three-way verdict (see interpret_rscape.py):"
echo "  significant pairs -> supports; none but adequate power -> rejects; none and low power -> cannot infer."
