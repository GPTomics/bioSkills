#!/bin/bash
# Reference: medaka 2.2+, seqkit 2.5+ | Verify API if version differs
# Polish an ONT-only draft assembly with medaka. ONE pass, no Racon pre-step,
# model auto-detected from the basecaller annotation in the reads.
set -euo pipefail

READS=${1:?Usage: $0 <reads.fq.gz> <draft.fa> [output_dir] [threads] [--bacteria]}
DRAFT=${2:?draft assembly required}
OUTPUT_DIR=${3:-medaka_output}
THREADS=${4:-8}            # parallelizes alignment/batching; the net itself is GPU-bound
BACTERIA=${5:-}           # pass --bacteria for native bacterial isolates (modified DNA)

# medaka is ONT-only. Do not run this on PacBio HiFi/CLR (no PacBio models exist).
# The model is auto-detected from the reads; supply -m only if auto-detection fails,
# and never copy a model name from an old tutorial onto new-chemistry data.
medaka_consensus -i "$READS" -d "$DRAFT" -o "$OUTPUT_DIR" -t "$THREADS" $BACTERIA

[ -f "${OUTPUT_DIR}/consensus.fasta" ] || { echo 'Error: polishing failed'; exit 1; }

echo 'Polishing complete. Draft vs polished length stats:'
seqkit stats "$DRAFT" "${OUTPUT_DIR}/consensus.fasta"
echo 'Validate the QV gain with reference-free Merqury on HELD-OUT k-mers,'
echo 'not the reads polished with (see genome-assembly/assembly-polishing).'
