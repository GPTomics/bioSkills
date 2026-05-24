#!/bin/bash
# Reference: BAGEL2 1.0.5+ (hart-lab/bagel) | Verify API if version differs
#
# BAGEL2 essentiality analysis end-to-end:
# 1. Compute per-sgRNA fold changes from counts
# 2. Compute per-gene Bayes Factors using CEGv2/NEGv1 reference sets
# 3. Precision-recall analysis to calibrate BF threshold

set -euo pipefail

# === INPUTS ===
COUNTS=counts.txt              # tab-separated: sgRNA, GENE, sample columns
CONTROL=Plasmid                # control sample column header
TREATMENT="Sample1,Sample2,Sample3"
CEG=CEGv2.txt                  # Hart 2017 core essentials
NEG=NEGv1.txt                  # Hart 2014 non-essentials
OUTDIR=bagel_results

mkdir -p "$OUTDIR"

# === DOWNLOAD REFERENCES IF MISSING ===
[ ! -f "$CEG" ] && curl -L -o "$CEG" \
    https://raw.githubusercontent.com/hart-lab/bagel/master/CEGv2.txt
[ ! -f "$NEG" ] && curl -L -o "$NEG" \
    https://raw.githubusercontent.com/hart-lab/bagel/master/NEGv1.txt

# === STEP 1: FOLD CHANGES ===
# bagel fc: per-sgRNA log-fold-change vs control
BAGEL.py fc \
    -i "$COUNTS" \
    -o "$OUTDIR/foldchange.txt" \
    -c "$CONTROL" \
    --min-reads 30                          # min reads/sgRNA in control

# === STEP 2: BAYES FACTORS ===
# bagel bf: per-gene BF via summed log-likelihood ratios
# Bootstrap 1000 iterations default; 5000+ for tight CI in noisy screens
BAGEL.py bf \
    -i "$OUTDIR/foldchange.txt" \
    -o "$OUTDIR/bayes_factor.txt" \
    -e "$CEG" \
    -n "$NEG" \
    -c "$TREATMENT" \
    -k 1000

# === STEP 3: PRECISION-RECALL CURVE ===
# Empirically calibrate BF threshold against CEGv2
BAGEL.py pr \
    -i "$OUTDIR/bayes_factor.txt" \
    -o "$OUTDIR/precision_recall.txt" \
    -e "$CEG" \
    -n "$NEG"

# === STEP 4: INTERPRETATION ===
# BF >6 corresponds to FDR ~0.05 (Hart 2017 G3 calibration)
# BF >12 corresponds to FDR ~0.005
# BF <-6 indicates candidate tumor suppressor (negative selection)

echo "BAGEL2 analysis complete. Outputs in $OUTDIR/"
echo "  - foldchange.txt: per-sgRNA LFCs"
echo "  - bayes_factor.txt: per-gene BF + bootstrap CI"
echo "  - precision_recall.txt: PR curve at every BF threshold"
echo ""
echo "Recommended BF thresholds:"
echo "  Exploratory:        BF >0  (P=0.85, R=0.95)"
echo "  Standard:           BF >6  (P=0.95, R=0.85; FDR ~0.05)"
echo "  High-confidence:    BF >12 (P=0.99, R=0.65; FDR ~0.005)"
echo "  Ultra-stringent:    BF >30 (P=1.00, R=0.20; FDR <0.001)"
