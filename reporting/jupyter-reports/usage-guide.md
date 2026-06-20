# Jupyter Reports Usage Guide

## Overview

This guide covers running parameterized Jupyter notebooks as batch report generators with papermill. The central idea is that a notebook carries hidden, mutable kernel state, so its saved outputs are not by themselves reproducible - the trustworthy artifact is the notebook re-executed top-to-bottom in a fresh kernel. papermill automates that clean-kernel execution and injects parameters; nbconvert renders the result to a code-hidden report. The decisions that matter are parameter handling, execution controls, how to aggregate results across samples, and pairing execution with a pinned environment so the report is reproducible in practice.

## Prerequisites

```bash
pip install papermill nbconvert scrapbook jupytext nbstripout

# Register the analysis environment as a kernel so papermill can target it
python -m ipykernel install --user --name analysis-env
```

## Quick Start

Tell your AI agent what you want to do:
- "Parameterize my DE-analysis notebook and run it on all 12 samples"
- "Generate a code-hidden HTML report from my analysis notebook"
- "Aggregate the per-sample QC metrics from my notebook runs into one table"
- "Stop my notebooks from bloating git with their outputs"
- "Make my notebook reproduce the same numbers on the cluster"

## Example Prompts

### Parameterizing and Batch Execution

> "Turn my QC notebook into a template that accepts sample_id and input_dir, and run it across the sample sheet"

> "Run my analysis template with three FDR thresholds and keep each executed notebook"

### Rendering

> "Convert my executed notebooks to HTML with the code hidden"

> "Make a PDF report from my notebook without installing LaTeX"

### Aggregation and Reproducibility

> "Collect the AUC and DE-gene count each notebook recorded into a cohort summary table"

> "Set up nbstripout so notebook outputs do not get committed, and pin the environment so reruns match"

## What the Agent Will Do

1. Add or confirm a single `parameters`-tagged cell with sensible defaults
2. Execute the template per sample in a fresh kernel with an explicit `execution_timeout`, preserving failed notebooks as artifacts
3. Render code-hidden HTML/PDF reports (`--no-input`, `--to webpdf` when LaTeX is unavailable)
4. Aggregate `scrapbook`-glued metrics across runs into one table
5. Wire nbstripout/jupytext for clean version control and pair execution with a pinned environment and seeded RNGs

## Tips

- Saved outputs are not reproducible; always re-execute (papermill or `nbconvert --execute`) rather than trusting them.
- Use `-r name value` (raw string) for parameters that look numeric but are identifiers (sample `007`, chromosome `1`); `-p` would coerce them.
- Set `execution_timeout` - the default is forever, which hangs a pipeline on one slow cell.
- The kernel must point at the pinned analysis environment; a registered-kernel mismatch is a top failure cause.
- papermill guarantees execution order, not the environment - pin packages (lockfile/container) and seed stochastic steps for true reproducibility.
- The notebook is the report, not the pipeline; move heavy multi-hour compute to Snakemake/Nextflow and keep the notebook for the summary.

## Related Skills

- reporting/quarto-reports - Document-first reporting that can render notebooks
- reporting/rmarkdown-reports - R-based literate reports
- workflows/scrnaseq-pipeline - Heavy compute as a pipeline, notebook as the summary
