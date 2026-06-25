#!/bin/bash
# Reference: Infernal 1.1.4+ | Verify API if version differs
# Build, calibrate, and search a custom covariance model from a structure-annotated alignment.

ALIGNMENT=$1
TARGET=$2
CM_NAME=${3:-"custom_family"}
THREADS=${4:-8}

if [ -z "$ALIGNMENT" ] || [ -z "$TARGET" ]; then
    echo "Usage: $0 <alignment.sto> <target.fa> [cm_name] [threads]"
    echo ""
    echo "The Stockholm alignment MUST include a #=GC SS_cons line; a structure-free alignment"
    echo "yields only an HMM-equivalent model. Validate the SS_cons covariation with R-scape"
    echo "(covariation-analysis) BEFORE building -- a non-covarying structure makes a confident-but-wrong CM."
    exit 1
fi

echo "=== Step 1: Build CM from the structure-annotated alignment ==="
cmbuild -n "$CM_NAME" "${CM_NAME}.cm" "$ALIGNMENT"

echo ""
echo "=== Step 2: Calibrate (required before E-values are meaningful; re-run after every rebuild) ==="
cmcalibrate --cpu "$THREADS" "${CM_NAME}.cm"

echo ""
echo "=== Step 3: Index ==="
cmpress "${CM_NAME}.cm"

echo ""
echo "=== Step 4: Search the target ==="
# E-value is valid only because Step 2 calibrated the model; without calibration switch to -T <bits>.
cmsearch \
    --cpu "$THREADS" \
    --tblout "${CM_NAME}_hits.tbl" \
    -E 1e-3 \
    "${CM_NAME}.cm" \
    "$TARGET" > "${CM_NAME}_hits.out"

echo ""
echo "=== Results ==="
NHITS=$(grep -cv '^#' "${CM_NAME}_hits.tbl" 2>/dev/null || echo 0)
echo "Hits found: $NHITS"
if [ "$NHITS" -gt 0 ]; then
    # Default cmsearch --tblout is fmt 1: score is column 15, E-value column 16.
    grep -v '^#' "${CM_NAME}_hits.tbl" | sort -k15,15 -rn | head -10 | \
        awk '{printf "  %s:%s-%s (%s) score=%.1f E=%s\n", $1, $8, $9, $10, $15, $16}'
fi

echo ""
echo "=== Optional: iterative refinement (watch for homology overextension) ==="
echo "  cmsearch -A new_hits.sto ${CM_NAME}.cm ${TARGET}"
echo "  # Curate new_hits.sto, merge with the seed, rebuild, then RE-CALIBRATE."
