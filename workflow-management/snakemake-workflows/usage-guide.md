# Snakemake Workflows - Usage Guide

## Overview

Snakemake is a pull, make-like workflow engine: an author declares target output files, each rule is a pattern-matched recipe, and the engine resolves the dependency DAG BACKWARD from the requested targets at parse time, before any job runs. The value of adopting it over a bash script is not "running steps in order" - it is reproducibility, provenance, and resumable caching. The engine buys reproducible workflow LOGIC automatically; the software environment, reference data, and hardware leaks are separate layers the author must pin (containers by digest, conda by lockfile). A clean DAG over unpinned tools is not reproducible. Reruns are decided by comparing on-disk state to a set of provenance triggers (mtime, params, input, code, software-env), and data-dependent branching - outputs whose number is unknown until a step runs - requires checkpoints, not imperative logic.

## Prerequisites

```bash
pip install snakemake

# Cluster execution (Snakemake 8+ executor plugin)
pip install snakemake-executor-plugin-slurm

# Cloud storage I/O (Snakemake 8+ storage plugin)
pip install snakemake-storage-plugin-s3

# Software deployment: conda (mamba recommended) and/or Apptainer/Singularity on PATH
```

Inputs: a `Snakefile` (or `workflow/Snakefile`), a `config/config.yaml`, per-rule `envs/*.yaml` or container images, and the analysis data. Run `snakemake --version` first - execution flags differ between Snakemake 7 and 8/9.

## Quick Start

Tell your AI agent what you want to do:
- "Create a Snakemake workflow for FASTQ-to-BAM alignment across all my samples"
- "Add a checkpoint so the pipeline scatters over however many contigs the assembler emits"
- "Port my Snakemake 7 `--cluster sbatch` command to the version 8 SLURM executor plugin"
- "Make each rule reproducible with a pinned conda env or container"
- "Escalate memory on retry so OOM-killed jobs get more RAM automatically"

## Example Prompts

### Authoring a Pipeline
> "Create a Snakemake workflow for RNA-seq: fastp trimming, salmon quantification, MultiQC, and a DESeq2 script, with wildcards fanning out over all samples in my config."

> "Wire a variant-calling pipeline where rules connect by output-file pattern, and add wildcard constraints so sample names containing underscores do not mis-route."

### Dynamic Outputs
> "The number of output files is unknown until a demultiplexing step runs - set up a checkpoint that splits into an unknown set of barcodes and an aggregation rule that gathers whatever was produced."

### Execution and HPC
> "Configure my workflow to run on SLURM with the version 8 executor plugin, mapping mem_mb and runtime to the scheduler, with a versioned profile."

> "My jobs are OOM-killed - add retries with memory that escalates on each attempt."

### Reproducibility and Debugging
> "Add pinned conda environments and digest-pinned containers to every rule and tell me how to run with software deployment enabled."

> "Explain why my whole pipeline re-ran after I edited a comment, and how to opt out of the code rerun trigger."

## What the Agent Will Do

1. Create a Snakefile with a `rule all` listing target outputs, and per-step rules whose output PATTERNS wire the DAG backward from those targets.
2. Define sample fan-out with wildcards + `expand()`, adding `wildcard_constraints` where values could contain delimiters or path separators.
3. Add checkpoints where the set of outputs is unknown until a step runs, with an input function that calls `checkpoints.X.get()` before globbing.
4. Declare `threads`, `resources` (mem_mb, runtime in minutes), and `retries` with attempt-based memory escalation for messy HPC jobs.
5. Pin software with per-rule `conda:` files and/or digest-pinned `container:` images, and select the deployment method (`--sdm conda apptainer`).
6. Branch execution guidance on the installed major version (7 vs 8/9) and provide dry-run, DAG, report, and cluster/cloud invocations.

## Key Concepts

| Concept | Description |
|---------|-------------|
| Rule | A pattern-matched recipe: an output pattern produced from inputs |
| Wildcard | A greedy regex token in a file pattern; constrain it to avoid mis-routing |
| expand() | Builds a list of target strings; does not touch the filesystem |
| Checkpoint | The only mechanism for outputs unknown until runtime (dynamic DAG) |
| Rerun trigger | mtime + params + input + code + software-env decide what recomputes |
| DAG | The static dependency graph built backward from requested targets |

## Tips

- Always `snakemake -n` (dry run) first - because the DAG is static, the dry run shows exactly what will run before any allocation is spent.
- Rules connect by OUTPUT-FILE PATTERN, not order in the file - a typo'd output path silently drops a rule; trace backward from `rule all` when a step does not run.
- Constrain wildcards (`wildcard_constraints: sample='[^/]+'`) whenever a value can contain `/`, `.`, or `_`; unconstrained greedy matching mis-routes silently.
- Checkpoints only re-evaluate the DAG through the exception raised by `checkpoints.X.get()` - glob the checkpoint's `directory()` output AFTER that call, never at parse time.
- Keep heavy compute out of `run:` blocks - they share the main process and GIL and cannot be isolated; use `script:` or `shell:`.
- Escalate memory on retry with `retries` + `mem_mb=lambda wildcards, attempt: base*attempt` rather than guessing one static value for the largest sample.
- On Snakemake 8/9 use `--executor slurm` (plugin) and `--sdm`, not the removed `--cluster`/`--use-conda`; check `snakemake --version` and branch.
- A pipeline is not reproducible until the environment is pinned - `conda:` files (or lockfiles) and containers pinned by `@sha256:` digest, never `:latest`.
- Use `--rerun-triggers mtime` to restore classic Make behavior when provenance triggers cause surprise reruns (or slow hashing on multi-TB inputs).
- For thousands of tiny HPC jobs, add `group:` so they submit as one job and scheduler latency stops dominating.

## Related Skills

- workflow-management/nextflow-pipelines - The reactive-dataflow alternative for runtime-emergent DAGs
- workflow-management/nf-core-pipelines - Run a curated community pipeline instead of authoring one
- workflows/rnaseq-to-de - An end-to-end RNA-seq pipeline this engine orchestrates
- read-qc/quality-reports - The QC step a pipeline wraps as an early rule
- read-alignment/bwa-alignment - The aligner invoked inside a Snakemake shell rule
