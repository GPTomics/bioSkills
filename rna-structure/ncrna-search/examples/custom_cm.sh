#!/bin/bash
# Build a custom covariance model from a Stockholm alignment and search a target database.

ALIGNMENT=$1
TARGET=$2
CM_NAME=${3:-"custom_family"}
THREADS=${4:-8}

if [ -z "$ALIGNMENT" ] || [ -z "$TARGET" ]; then
    echo "Usage: $0 <alignment.sto> <target.fa> [cm_name] [threads]"
    echo ""
    echo "Builds a CM from a structure-annotated Stockholm alignment"
    echo "and searches a target database."
    echo ""
    echo "The Stockholm alignment must include #=GC SS_cons line."
    exit 1
fi

echo "=== Step 1: Build CM from alignment ==="
cmbuild -n "$CM_NAME" "${CM_NAME}.cm" "$ALIGNMENT"

echo ""
echo "=== Step 2: Calibrate CM (this may take several minutes) ==="
cmcalibrate --cpu "$THREADS" "${CM_NAME}.cm"

echo ""
echo "=== Step 3: Index CM ==="
cmpress "${CM_NAME}.cm"

echo ""
echo "=== Step 4: Search target database ==="
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
    echo ""
    echo "Top hits by score:"
    grep -v '^#' "${CM_NAME}_hits.tbl" | sort -k15,15 -rn | head -10 | \
        awk '{printf "  %s:%s-%s (%s) score=%.1f E=%s\n", $1, $8, $9, $10, $15, $16}'
fi

echo ""
echo "=== Optional: Iterative refinement ==="
echo "To refine the CM with new hits:"
echo "  cmsearch -A new_hits.sto ${CM_NAME}.cm ${TARGET}"
echo "  # Manually curate new_hits.sto, merge with original alignment"
echo "  # Rebuild: cmbuild -n ${CM_NAME} ${CM_NAME}_v2.cm merged.sto"
