# Automated QC Reports Usage Guide

## Overview

This guide covers aggregating per-tool QC metrics into one MultiQC report. The central idea: MultiQC does not measure anything - it scrapes the metrics other tools already wrote and presents them. The report is a human-triage snapshot, not a pass/fail gate, so the real decisions are how to scope detection, how to keep sample names from colliding, and how to turn the parsed data into an actual gate when one is needed.

## Prerequisites

```bash
pip install multiqc

# Or via conda
conda install -c bioconda multiqc
```

Generate the upstream QC outputs first (FastQC/falco, fastp, STAR/Picard logs, salmon, bcftools, etc.). MultiQC has nothing to aggregate without them.

## Quick Start

Tell your AI agent what you want to do:
- "Generate a MultiQC report from my FastQC and alignment results"
- "Aggregate QC metrics across all samples and tell me which ones look off"
- "My MultiQC report merged two samples into one row - fix the sample names"
- "Build a QC gate that fails the pipeline when mapping rate drops below 70%"
- "Pin a reproducible MultiQC config for my pipeline"

## Example Prompts

### Basic Reports

> "Run MultiQC on my results directory and produce a shareable HTML report"

> "Combine FastQC, STAR, and featureCounts outputs into one QC summary"

### Scoping and Sample Names

> "MultiQC is picking up files that aren't real tool outputs - scope it to just the modules I ran"

> "Two lanes of the same sample collapsed into one row; configure sample-name cleaning so they stay separate"

### From Report to Gate

> "Parse multiqc_data.json and fail the job if any sample has under 50% mapped reads"

> "Set up MultiQC in my Snakemake workflow with a pinned config and a separate QC-gating step"

### Governance

> "Generate the report but make sure no sample names or data leave our network"

## What the Agent Will Do

1. Confirm which upstream QC outputs exist and check the expected sample roster
2. Scope detection with `-m`/`--ignore`/`sp:` so only intended modules and files are parsed
3. Resolve sample-name collisions via `extra_fn_clean_exts`, `--fn_as_s_name`, `--dirs`, or `--replace-names`, verifying against `multiqc_sources.txt`
4. Pin a `multiqc_config.yaml` (title, `module_order`, formatting) for reproducible reports
5. When gating is requested, emit `multiqc_data.json` and apply thresholds in a separate parse-and-exit step rather than expecting MultiQC to fail on bad QC
6. Disable AI summaries for sensitive data

## Tips

- An empty report usually means "nothing matched," not "QC passed" - check the header sample count against what was expected.
- Read `multiqc_data/multiqc_sources.txt` first when rows are merged, split, or missing; it maps every file to the sample name it produced.
- Conditional-formatting colors are configured thresholds, not biological verdicts - own the thresholds before trusting the colors.
- `--data-format json` gives the machine-readable `multiqc_data.json`; build automation on that, not on the HTML.
- Pin the MultiQC version alongside the tool versions - parsers and defaults shift between releases.
- AI summaries are off by default and send metrics plus sample names off-network when enabled; use `--no-ai` for clinical or embargoed data.

## Related Skills

- read-qc/quality-reports - Generate the FastQC/falco inputs MultiQC aggregates
- read-qc/fastp-workflow - Preprocessing QC that feeds the report
- workflows/rnaseq-to-de - Full workflow that emits a MultiQC report
