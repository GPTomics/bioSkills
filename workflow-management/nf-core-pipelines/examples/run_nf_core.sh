#!/usr/bin/env bash
# Reference: nf-core/Nextflow 24.04+ | Verify API if version differs
# Runs nf-core/rnaseq reproducibly: pinned -r release, comma-separated -profile,
# a validated samplesheet, an explicit outdir, and -resume. Verify flags with
# `nextflow run nf-core/rnaseq -r <ver> --help` if the installed version differs.
set -euo pipefail

PIPELINE='nf-core/rnaseq'
REVISION='3.14.0'                 # pin the release; a mutable default branch is not reproducible
PROFILE='singularity'            # container engine; comma-add an institution, e.g. singularity,uppmax
SAMPLESHEET='samplesheet.csv'    # columns must match the pipeline schema (assets/schema_input.json)
GENOME='GRCh38'                  # iGenomes key; swap for explicit --fasta/--gtf to pin a build
OUTDIR='results'

# Smoke-test first: a tiny public dataset proves the install + engine end to end.
nextflow run "$PIPELINE" -r "$REVISION" -profile "test,${PROFILE}" --outdir "${OUTDIR}_test"

# Real run: single-dash options are Nextflow, double-dash are pipeline parameters.
nextflow run "$PIPELINE" -r "$REVISION" \
    -profile "$PROFILE" \
    --input "$SAMPLESHEET" \
    --genome "$GENOME" \
    --outdir "$OUTDIR" \
    -resume

echo "Done. Triage per-sample QC in ${OUTDIR}/multiqc/star_salmon/multiqc_report.html before trusting results."
