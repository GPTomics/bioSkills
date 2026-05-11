#!/usr/bin/env bash
# Run nf-core/rnaseq on LANTA with Slurm + Apptainer.
# Launch from lanta.nstda.or.th after staging data/containers via transfer.lanta.nstda.or.th.

set -euo pipefail

export PROJECT_ID="${PROJECT_ID:-<PROJECT_ID-SHORTNAME>}"
export LANTA_ACCOUNT="${LANTA_ACCOUNT:-<LANTA_ACCOUNT>}"
NFCORE_RNASEQ_VERSION="${NFCORE_RNASEQ_VERSION:-3.25.0}"
SAMPLESHEET="${SAMPLESHEET:-rnaseq_samplesheet.csv}"
FASTA="${FASTA:-/project/${PROJECT_ID}/refs/genome.fa}"
GTF="${GTF:-/project/${PROJECT_ID}/refs/genes.gtf}"
export OUTDIR="${OUTDIR:-/project/${PROJECT_ID}/results/rnaseq}"
export WORKDIR="${WORKDIR:-/scratch/${PROJECT_ID}/${USER}/nextflow-work/rnaseq}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-/scratch/${PROJECT_ID}/${USER}/apptainer-cache}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-/scratch/${PROJECT_ID}/${USER}/apptainer-tmp}"
CONFIG="${CONFIG:-nextflow.lanta.config}"
NEXTFLOW_MODULE="${NEXTFLOW_MODULE:-}"
APPTAINER_MODULE="${APPTAINER_MODULE:-}"

if [[ "$PROJECT_ID" == "<PROJECT_ID-SHORTNAME>" || "$LANTA_ACCOUNT" == "<LANTA_ACCOUNT>" ]]; then
  echo "ERROR: Set PROJECT_ID and LANTA_ACCOUNT before running." >&2
  exit 1
fi

if [[ ! -f "$SAMPLESHEET" ]]; then
  echo "ERROR: samplesheet not found: $SAMPLESHEET" >&2
  exit 1
fi

if [[ ! -f "$FASTA" ]]; then
  echo "ERROR: FASTA not found: $FASTA" >&2
  exit 1
fi

if [[ ! -f "$GTF" ]]; then
  echo "ERROR: GTF not found: $GTF" >&2
  exit 1
fi

mkdir -p "$OUTDIR" "$WORKDIR" "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" logs

if command -v module >/dev/null 2>&1; then
  module purge || true
  if [[ -n "$NEXTFLOW_MODULE" ]]; then module load "$NEXTFLOW_MODULE"; fi
  if [[ -n "$APPTAINER_MODULE" ]]; then module load "$APPTAINER_MODULE"; fi
fi

command -v nextflow >/dev/null 2>&1 || {
  echo "ERROR: nextflow not found. On LANTA, check module avail nextflow and set NEXTFLOW_MODULE." >&2
  exit 1
}

if ! command -v apptainer >/dev/null 2>&1 && ! command -v singularity >/dev/null 2>&1; then
  echo "ERROR: Apptainer/Singularity not found. On LANTA, check module avail apptainer/singularity and set APPTAINER_MODULE." >&2
  exit 1
fi

nextflow run nf-core/rnaseq \
  -r "$NFCORE_RNASEQ_VERSION" \
  -c "$CONFIG" \
  --input "$SAMPLESHEET" \
  --outdir "$OUTDIR" \
  --fasta "$FASTA" \
  --gtf "$GTF" \
  -work-dir "$WORKDIR" \
  -resume \
  -with-report "$OUTDIR/report.html" \
  -with-trace "$OUTDIR/trace.txt" \
  -with-timeline "$OUTDIR/timeline.html" \
  -with-dag "$OUTDIR/flow.svg"

echo "Submitted/ran nf-core/rnaseq on LANTA. Check $OUTDIR and .nextflow.log. Use myqueue for Slurm status."
