# nf-core Pipelines - Usage Guide

## Overview

nf-core is a curated collection of community-built Nextflow pipelines (rnaseq, sarek, atacseq, methylseq, ampliseq, taxprofiler, fetchngs, and dozens more) that encode years of QC, edge-case handling, CI/nf-test regression tests, institutional configs, standardized samplesheets, and MultiQC reporting. The governing decision for most biologists is ADOPT a community pipeline rather than author one from scratch: reinventing a mainstream analysis in hand-written Nextflow is months of work and worse QC. The second thing to internalize is that adoption only delivers reproducibility if the RUN is pinned - `-r` pins the pipeline version, the release's digest-pinned containers pin the software, and explicit references pin the reference data. An unpinned run on the default branch with a `:latest` engine is exactly as irreproducible as a bespoke script, so the curation buys nothing unless the invocation is disciplined.

## Prerequisites

```bash
# Nextflow (JVM required: Java 17+)
curl -s https://get.nextflow.io | bash

# nf-core helper tools (list/download pipelines, scaffold/lint your own)
pip install nf-core

# A container engine (pick the one your platform allows)
#   Docker 24+          - laptop/workstation with root/daemon
#   Singularity/Apptainer 3.8+ - shared HPC, rootless (most common on clusters)
#   conda               - last resort when no container engine is available
```

Inputs: a samplesheet CSV whose columns match the target pipeline's schema, reference data (an iGenomes `--genome` key or explicit `--fasta`/`--gtf`), and optionally a custom or institutional config for HPC/cloud execution.

## Quick Start

Tell your AI agent what you want to do:
- "Run the nf-core/rnaseq test profile to confirm my install works"
- "Run nf-core/rnaseq on my paired-end samples pinned to release 3.14.0 with Singularity"
- "Build a samplesheet for nf-core/rnaseq from my FASTQ directory"
- "Run nf-core/sarek for somatic variant calling on tumor/normal pairs"
- "Configure nf-core/rnaseq to run on my SLURM cluster with an institutional config"
- "Read the MultiQC report and flag any samples that failed QC"

## Example Prompts

### Adopt and Run a Pipeline
> "Run nf-core/rnaseq pinned to a specific release on my paired-end human samples using Singularity, with GRCh38, and tell me at the end which samples look bad in MultiQC."

### Build a Samplesheet
> "I have a directory of paired FASTQ files named SAMPLE_R1.fastq.gz / SAMPLE_R2.fastq.gz. Build a valid nf-core/rnaseq samplesheet.csv with the correct columns and validate it against the pipeline schema before running."

### Somatic Variant Calling
> "Run nf-core/sarek for somatic calling on three tumor/normal pairs. Set up the patient/sample/lane samplesheet, pick the right variant callers, and pin the pipeline revision."

### HPC / Institutional Config
> "Configure and launch nf-core/atacseq on my SLURM cluster. Use an institutional config from nf-core/configs if one exists for my site, otherwise write a custom config that sets the executor, queue, and max_memory, and enable resume."

### Reference Choice
> "Should I use --genome GRCh38 or supply my own --fasta and --gtf for this run? I need a specific Ensembl release - set it up the reproducible way and record the reference source."

## What the Agent Will Do

1. Confirm ADOPT-vs-BUILD: check whether a community pipeline covers the analysis, and default to adopting it rather than authoring Nextflow.
2. Select the pipeline and pin a specific `-r` release for reproducibility.
3. Choose a `-profile`: a container engine (docker/singularity/conda) plus, where available, an institutional config from nf-core/configs, comma-separated with no spaces.
4. Build the samplesheet CSV with the pipeline's exact schema columns and validate it (a `-stub` or `-profile test` pre-flight) before committing compute.
5. Set references via `--genome` (iGenomes) or explicit `--fasta`/`--gtf`, recording the source for provenance.
6. Configure resources for the execution platform, capping retry escalation with `max_memory`/`max_cpus`/`max_time`.
7. Launch with `-resume`, then read the aggregated `multiqc_report.html` to triage per-sample QC before trusting downstream results.

## Tips

- Always pin `-r <version>`; without it the run tracks the mutable default branch and cannot be reproduced.
- Profiles are comma-separated with NO spaces: `-profile test,docker`, not `-profile test docker`.
- Run `-profile test,<engine>` (a tiny public dataset) or `-stub` once before real data to prove the install end to end.
- Pipeline parameters take a double dash (`--input`, `--genome`, `--outdir`); Nextflow options take a single dash (`-r`, `-profile`, `-resume`) - mixing them up is the most common error.
- Match the samplesheet columns to the pipeline's own schema (`docs/usage.md` / `assets/schema_input.json`); the columns differ per pipeline (rnaseq vs sarek).
- On HPC, enable `-resume` and add `cache 'lenient'` via a custom config so network-filesystem mtime quirks do not force full re-runs; never delete `work/` if you might resume.
- Prefer explicit `--fasta`/`--gtf` over `--genome` when you need a specific build - iGenomes references are frozen and can lag current annotation releases.
- Layer overrides through `-c custom.config`; do not edit the pipeline's `conf/base.config` so the checkout stays clean and updatable.
- Read MultiQC as a triage surface, not a verdict: it reports what tools measured; you set the pass/fail thresholds from the assay.

## Related Skills

- workflow-management/nextflow-pipelines - Author your own Nextflow pipeline when no community pipeline fits
- workflow-management/snakemake-workflows - Rule-based alternative engine for pipeline authoring
- workflows/rnaseq-to-de - Take an nf-core/rnaseq count matrix into differential expression
- read-qc/quality-reports - Interpret the FastQC/MultiQC QC surface a pipeline emits
