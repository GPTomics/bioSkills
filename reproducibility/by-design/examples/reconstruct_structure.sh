#!/usr/bin/env bash
# Reference: BioFlow-Insight 0.2+, Nextflow 24.04+, Snakemake 8+, Graphviz 2.50+ | Verify CLI if version differs
#
# Reproducibility-by-design, pillar 1 (modularity) + comprehension: derive a workflow's
# structure diagram FROM THE CODE so it can never drift from what actually runs. A hand-drawn
# slide looks the same but decays; a code-derived graph is regenerated on every change.
#
# Usage: ./reconstruct_structure.sh <path-to-pipeline>
#   where <path-to-pipeline> is a Nextflow main.nf OR a Snakemake Snakefile.
set -euo pipefail

PIPELINE="${1:?usage: reconstruct_structure.sh <main.nf | Snakefile>}"
OUTDIR="structure"
mkdir -p "$OUTDIR"

case "$PIPELINE" in
  *.nf)
    echo "[nextflow] reconstructing structure from Nextflow DSL2 code"

    # 1) BioFlow-Insight: parse the AST and reconstruct the workflow structure (Marchment 2024).
    #    The graph is DERIVED from the code, so clicking a step maps back to the exact lines.
    if command -v bioflow-insight >/dev/null 2>&1; then
      bioflow-insight "$PIPELINE" --output-dir "$OUTDIR/" || \
        echo "  (confirm current flags with: bioflow-insight --help)"
    else
      echo "  bioflow-insight not installed: pip install bioflow-insight"
    fi

    # 2) Nextflow's own code-derived DAG. -preview parses WITHOUT running any task,
    #    so this is safe on any machine and needs no input data.
    if command -v nextflow >/dev/null 2>&1; then
      nextflow run "$PIPELINE" -preview -with-dag "$OUTDIR/flowchart.mmd" || true
      echo "  wrote $OUTDIR/flowchart.mmd (Mermaid; also .dot/.svg/.png by extension)"
    fi
    ;;

  *)
    echo "[snakemake] reconstructing structure from Snakemake code"
    # --rulegraph is the module-level view (one node per rule), closest to reusable steps.
    # --dag is the per-job expansion. Both are generated from the code, never hand-drawn.
    snakemake -s "$PIPELINE" --rulegraph | dot -Tsvg > "$OUTDIR/rulegraph.svg"
    snakemake -s "$PIPELINE" --dag       | dot -Tsvg > "$OUTDIR/dag.svg"
    echo "  wrote $OUTDIR/rulegraph.svg (modules) and $OUTDIR/dag.svg (per-job)"
    ;;
esac

echo
echo "Design check (pillar 1): for each step in the diagram, can a consumer answer"
echo "  (a) what it does, (b) how to use it, (c) what it depends on"
echo "without reading the internals? If not, the module boundary needs work."
echo "Commit the regenerated diagram alongside the code so it never drifts."
