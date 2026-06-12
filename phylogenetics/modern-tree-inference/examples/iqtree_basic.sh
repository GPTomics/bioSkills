#!/bin/bash
# Reference: IQ-TREE 2.2+/2.3+ | Verify API if version differs
# Standard ML tree: ModelFinder + dual support (UFBoot2 + SH-aLRT).
# Reframe: the support is repeatability under resampling, NOT correctness.
# A branch is strongly supported only if SH-aLRT >= 80 AND UFBoot >= 95.
# NOT spot-runnable offline: needs the iqtree2 binary and a real alignment.
set -euo pipefail

ALIGNMENT="${1:-alignment.fasta}"
OUT="iqtree_out"                       # all outputs namespaced here, never the CWD
mkdir -p "$OUT"

# -m MFP   ModelFinder Plus: BIC-rank all models (incl. FreeRate +R), then search
# -B 1000  ultrafast bootstrap, >=1000 reps          (v1.x used -bb)
# -bnni    UFBoot2 NNI safeguard against model-violation inflation; pair with -B
# -alrt 1000  SH-aLRT, the tree-perturbation companion to UFBoot
# -T AUTO  auto thread count; -ntmax caps it          (v1.x used -nt)
# --seed   reproducibility
iqtree2 -s "$ALIGNMENT" -m MFP -B 1000 -bnni -alrt 1000 -T AUTO -ntmax 8 \
        --seed 12345 --prefix "$OUT/run1"

echo "Best tree:  $OUT/run1.treefile"
echo "Report:     $OUT/run1.iqtree"
grep "Best-fit model" "$OUT/run1.iqtree" || true

# Node labels are written as SH-aLRT/UFBoot (e.g. 92.5/98). Strong iff both pass:
echo "Interpret labels as SH-aLRT/UFBoot; strong support needs >=80 AND >=95."
