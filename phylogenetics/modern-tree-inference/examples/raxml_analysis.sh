#!/bin/bash
# Reference: RAxML-NG 1.2+ | Verify API if version differs
# RAxML-NG for very large trees: ML search + bootstrap with BOTH the Felsenstein
# proportion and the transfer bootstrap (TBE), which rescues deep branches that a
# single rogue taxon would crush in a big tree.
# NOT spot-runnable offline: needs the raxml-ng binary and a real alignment.
set -euo pipefail

ALIGNMENT="${1:-alignment.fasta}"
OUT="raxml_out"
mkdir -p "$OUT"

# Pre-flight: validate the alignment and estimate RAM/threads before the run
raxml-ng --check --msa "$ALIGNMENT" --model GTR+G --prefix "$OUT/check"
raxml-ng --parse --msa "$ALIGNMENT" --model GTR+G --prefix "$OUT/parse"

# --all                ML search + bootstrap + draw support in one command
# --bs-trees autoMRE{1000}  bootstrap with MRE convergence test, cap 1000
# --bs-metric fbp,tbe  Felsenstein proportion AND transfer bootstrap expectation
# --threads auto{8}    auto thread/worker detection capped at 8
raxml-ng --all --msa "$ALIGNMENT" --model GTR+G \
         --bs-trees 'autoMRE{1000}' --bs-metric fbp,tbe \
         --threads 'auto{8}' --seed 42 --prefix "$OUT/run_ng"

echo "Best tree:    $OUT/run_ng.raxml.bestTree"
echo "With support: $OUT/run_ng.raxml.support (FBP and TBE columns)"

# Separate steps for a long HPC run with checkpointing:
# raxml-ng --search    --msa "$ALIGNMENT" --model GTR+G --seed 42 --prefix "$OUT/ml"
# raxml-ng --bootstrap --msa "$ALIGNMENT" --model GTR+G --seed 42 --bs-trees 1000 --prefix "$OUT/bs"
# raxml-ng --support   --tree "$OUT/ml.raxml.bestTree" --bs-trees "$OUT/bs.raxml.bootstraps" \
#          --bs-metric tbe --prefix "$OUT/sup"
