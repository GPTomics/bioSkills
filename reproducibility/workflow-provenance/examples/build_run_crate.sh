#!/usr/bin/env bash
# Reference: cwltool 3.1+, runcrate 0.5+, rocrate 0.10+, Nextflow 24.04+, nf-prov 1.2+ | Verify CLI if version differs
#
# Capture the provenance of ONE workflow run and package it as a portable, citable
# Workflow Run RO-Crate. Two reproducibility-by-design points are enforced:
#   1) provenance != reproducibility: the RO records lineage, but only pinning the
#      environment (digests) makes the recorded run reproducible;
#   2) use the Provenance/Run profile (records the execution), not the plain Workflow profile.
#
# Usage: ./build_run_crate.sh cwl   workflow.cwl inputs.yml
#        ./build_run_crate.sh nextflow main.nf   (nf-prov must be enabled in nextflow.config)
set -euo pipefail

ENGINE="${1:?usage: build_run_crate.sh <cwl|nextflow> ...}"; shift

case "$ENGINE" in
  cwl)
    WF="${1:?workflow.cwl}"; JOB="${2:?inputs.yml}"
    # 1) Run with provenance ON -> a CWLProv Research Object (PROV + inputs/outputs/intermediates).
    cwltool --provenance prov-ro/ "$WF" "$JOB"
    # 2) Convert the CWLProv RO into a Workflow Run RO-Crate and print the lineage.
    runcrate convert prov-ro/ --output run-crate/
    runcrate report  run-crate/
    echo "  Workflow Run RO-Crate at run-crate/ (ro-crate-metadata.json is JSON-LD PROV)"
    ;;

  nextflow)
    WF="${1:?main.nf}"
    # nextflow.config must enable: plugins { id 'nf-prov' } and a trace/report/dag (see SKILL.md).
    # nf-prov emits a Workflow Run RO-Crate (and can emit a BioCompute Object) for the run.
    nextflow run "$WF" -with-trace prov/trace.txt -with-report prov/report.html -with-dag prov/dag.mmd
    echo "  trace/report/dag under prov/ ; nf-prov RO-Crate at prov/ro-crate-metadata.json"
    echo "  REPRODUCIBILITY CHECK: grep the trace for ':latest' -- any moving tag means the"
    echo "  recorded run is NOT reproducible. Pin every container by @sha256 digest and re-run."
    grep -q ':latest' prov/trace.txt 2>/dev/null && echo "  WARNING: ':latest' found in trace" || true
    ;;

  *)
    echo "unknown engine: $ENGINE (use cwl or nextflow)"; exit 1 ;;
esac

echo
echo "Publish: register run-crate/ (or prov/) on WorkflowHub for a citable DOI,"
echo "or on Dockstore (GA4GH TRS) to make it discoverable and runnable across platforms."
