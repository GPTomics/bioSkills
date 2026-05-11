#!/usr/bin/env bash
# Run nf-core/rnaseq on RAPTOR / MedCMU-HPC with Slurm + Apptainer.
# Launch from the RAPTOR login node; compute work is submitted through Slurm by Nextflow.

set -euo pipefail

RAPTOR_PARTITION="${RAPTOR_PARTITION:-<PARTITION>}"
PROJECT_DIR="${PROJECT_DIR:-<PROJECT_DIR>}"
SCRATCH_DIR="${SCRATCH_DIR:-<SCRATCH_DIR>}"
NFCORE_RNASEQ_VERSION="${NFCORE_RNASEQ_VERSION:-3.25.0}"
SAMPLESHEET="${SAMPLESHEET:-rnaseq_samplesheet.csv}"
FASTA="${FASTA:-${PROJECT_DIR}/refs/genome.fa}"
GTF="${GTF:-${PROJECT_DIR}/refs/genes.gtf}"
OUTDIR="${OUTDIR:-${PROJECT_DIR}/results/rnaseq}"
WORKDIR="${WORKDIR:-${SCRATCH_DIR}/nextflow-work/rnaseq}"
APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${SCRATCH_DIR}/apptainer-cache}"
APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-${SCRATCH_DIR}/apptainer-tmp}"
CONFIG="${CONFIG:-nextflow.raptor.config}"
NEXTFLOW_MODULE="${NEXTFLOW_MODULE:-nextflow}"
APPTAINER_MODULE="${APPTAINER_MODULE:-apptainer}"

if [[ "$RAPTOR_PARTITION" == "<PARTITION>" || "$PROJECT_DIR" == "<PROJECT_DIR>" || "$SCRATCH_DIR" == "<SCRATCH_DIR>" ]]; then
  echo "ERROR: Set RAPTOR_PARTITION, PROJECT_DIR, and SCRATCH_DIR before running." >&2
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
  module load "$APPTAINER_MODULE"
  module load "$NEXTFLOW_MODULE"
fi

command -v nextflow >/dev/null 2>&1 || {
  echo "ERROR: nextflow not found. Check module avail nextflow on RAPTOR." >&2
  exit 1
}

command -v apptainer >/dev/null 2>&1 || {
  echo "ERROR: apptainer not found. Check module avail apptainer on RAPTOR." >&2
  exit 1
}

export APPTAINER_CACHEDIR APPTAINER_TMPDIR

nextflow run nf-core/rnaseq \
  -r "$NFCORE_RNASEQ_VERSION" \
  -c "$CONFIG" \
  --input "$SAMPLESHEET" \
  --outdir "$OUTDIR" \
  --fasta "$FASTA" \
  --gtf "$GTF" \
  --raptor_partition "$RAPTOR_PARTITION" \
  --project_dir "$PROJECT_DIR" \
  --scratch_dir "$SCRATCH_DIR" \
  --apptainer_cache "$APPTAINER_CACHEDIR" \
  --apptainer_tmp "$APPTAINER_TMPDIR" \
  -work-dir "$WORKDIR" \
  -resume \
  -with-report "$OUTDIR/report.html" \
  -with-trace "$OUTDIR/trace.txt" \
  -with-timeline "$OUTDIR/timeline.html" \
  -with-dag "$OUTDIR/flow.svg"

echo "Submitted/ran nf-core/rnaseq on RAPTOR. Check $OUTDIR, .nextflow.log, squeue, and sacct."
