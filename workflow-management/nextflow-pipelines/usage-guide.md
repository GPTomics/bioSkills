# Nextflow Pipelines - Usage Guide

## Overview

Nextflow is a reactive-dataflow workflow engine: processes are wired together by asynchronous channels, each channel item is a future, and a task fires the instant its inputs are ready, so execution order is not guaranteed. This is the opposite of pull/goal-oriented engines (Snakemake/CWL/WDL) that reason backward from a target file. The payoff is reproducibility, resumable caching, and executor portability (the same code moves from a laptop to SLURM to AWS Batch by swapping a profile). The trap to internalize first is that the DAG buys reproducibility of workflow LOGIC only; a pipeline over unpinned tools is not reproducible, so the software environment (container by digest, conda by lockfile), reference data, and thread/locale/arch leaks must be pinned by hand. Most of the field's day-to-day pain (only the first sample ran; `-resume` re-ran everything) traces directly back to the dataflow model and to unpinned or nondeterministic inputs.

## Prerequisites

```bash
# Install Nextflow (requires a Java 17+ JVM)
curl -s https://get.nextflow.io | bash

# Or via conda
conda install -c bioconda nextflow

# A container runtime for reproducible per-process software:
#   Docker (local/cloud) or Singularity/Apptainer (HPC, rootless)
```

Inputs: a `main.nf` (DSL2), an optional `nextflow.config` with profiles, per-process containers pinned by immutable tag or digest, and read/reference paths supplied via `params` or a samplesheet.

## Quick Start

Tell your AI agent what you want to do:
- "Author a Nextflow DSL2 pipeline for paired-end RNA-seq quantification"
- "Fix my pipeline so it processes every sample, not just the first"
- "Diagnose why my `-resume` keeps re-running everything on the cluster"
- "Add dynamic memory retry escalation to my alignment process"
- "Make my pipeline portable across local, SLURM, and AWS Batch with profiles"

## Example Prompts

### Authoring a Pipeline
> "Build a Nextflow DSL2 pipeline with FASTQC, fastp trimming, Salmon quantification, and MultiQC, using per-process containers pinned by tag and a meta-map convention for sample identity."

> "Refactor my monolithic main.nf into DSL2 modules and a QC subworkflow with take/main/emit."

### Channels and Dataflow
> "My pipeline only processes the first sample and exits cleanly with no error - explain why and fix it."

> "Set up channels to join per-sample BAMs with their metadata by sample key without silently dropping unmatched samples."

### Resume and Caching
> "My `-resume` works locally but misses the cache on our Lustre cluster and re-runs everything - diagnose and fix it."

> "Use `-dump-hashes` and `nextflow log` to find which task's hash changed between two runs."

### Resources and Portability
> "Add exit-status-conditional retry with memory escalation to my OOM-prone process."

> "Configure profiles so the same pipeline runs on local Docker, SLURM with Singularity, and AWS Batch."

## What the Agent Will Do

1. Author a `main.nf` in DSL2 with processes wired through channels, and split reusable tools into modules with take/main/emit subworkflows.
2. Choose queue vs value channels correctly, using `.first()` (or `Channel.value`) for any shared reference so every sample fires, not just the first.
3. Select operators with their traps in mind (`groupTuple` size hint, `join` remainder to avoid silent drops).
4. Pin per-process containers by immutable tag or digest and set the software runtime and executor in profiles, never in pipeline code.
5. Add dynamic retry escalation (`memory { 8.GB * task.attempt }`) with exit-status-conditional `errorStrategy`.
6. Diagnose `-resume` cache misses (nondeterministic ordering, mtime, mutable tags) with `cache 'lenient'`, `-dump-hashes`, and `nextflow log`.
7. Recommend adopting an nf-core community pipeline when the analysis is mainstream, reserving authoring for genuinely novel logic.

## Tips

- A shared reference must come from a VALUE channel (`.first()`/`Channel.value`); a queue channel is drained after the first task and later samples silently never run.
- `groupTuple` without `size:` or `groupKey` waits for the entire upstream channel to close - one slow sample stalls all grouping; a never-closing channel deadlocks.
- `join` silently drops non-matching keys by default; use `remainder: true` or `failOnMismatch: true` when a missing key means a bug, not a design choice.
- `work/<hash>/` is the pipeline's real output store and the only thing `-resume` reads; `publishDir` is a side-effect copy whose failure can be silent. Never `rm -rf work/`; use `nextflow clean`.
- Pin containers by immutable digest; `:latest` breaks reproducibility and can serve a stale cache hit (wrong result, no error) after the image moves.
- On network filesystems set `cache 'lenient'` so mtime differences do not bust the cache; this is the single most common HPC resume fix.
- Sort `collect`/`groupTuple` outputs before a downstream task, or their nondeterministic order changes the task hash and forces a re-run.
- Escalate resources per attempt (`memory { 8.GB * task.attempt }`) and gate retries on OOM/kill exit codes so transient failures auto-recover without wasting retries on deterministic errors.
- Keep the executor in a profile, never in the pipeline body, so the same code runs local, on SLURM/LSF, and on AWS/Google Batch.
- Adopt `nf-core/<pipeline>` with a pinned `-r` revision before authoring anything mainstream; DIY is a permanent maintenance liability.

## Related Skills

- workflow-management/nf-core-pipelines - Run and configure community pipelines (the running counterpart to this authoring skill)
- workflow-management/snakemake-workflows - Pull/goal-oriented alternative for Python shops and file-pattern logic
- workflows/rnaseq-to-de - End-to-end RNA-seq quantification to differential expression
- read-qc/quality-reports - QC steps a pipeline orchestrates (FastQC/MultiQC)
