# workflow-management

## Overview

Workflow engines for authoring bioinformatics pipelines. Their value is not "running steps in order" - a bash script does that - but reproducibility, provenance, and resumable caching: a declarative DAG the engine schedules, content/parameter-hash caching so a re-run recomputes only what changed, and a pinned software environment so the result is the same next year on another machine. The critical caveat: adopting an engine buys reproducible workflow LOGIC only; the author must still pin the software environment (containers by DIGEST not mutable tag, conda by LOCKFILE), the reference data and parameters, and control hardware/OS leaks (thread count, sort order, CPU architecture). A clean DAG over unpinned tools is not reproducible.

**Tool type:** mixed | **Primary tools:** Snakemake, Nextflow, nf-core, cwltool, Cromwell/miniwdl

## Skills

| Skill | Description |
|-------|-------------|
| snakemake-workflows | Author pull/goal-oriented pipelines with Snakemake rules, wildcards, checkpoints, and the v8 executor-plugin model |
| nextflow-pipelines | Author reactive-dataflow pipelines with Nextflow DSL2 channels, resume caching, and executor portability |
| nf-core-pipelines | Run and configure curated nf-core community pipelines with pinned revisions, samplesheets, and institutional configs |
| cwl-workflows | Author portable, standards-based pipelines with the Common Workflow Language spec and secondaryFiles |
| wdl-workflows | Author WDL pipelines for the Terra/AnVIL and GATK/Broad ecosystem with cost-aware runtime blocks |

## Choosing an engine

The choice is driven by the ecosystem to integrate with, not by which engine is "best". The first axis is the execution model: pull/goal-oriented engines (Snakemake, CWL, WDL) declare target outputs and build a mostly-static DAG backward, so a dry-run shows the plan before anything runs; Nextflow is push/reactive-dataflow, wiring processes through channels where the DAG emerges at runtime and dynamic branching is native.

| Engine | Model | Sweet spot | Community pipelines | When to choose it |
|--------|-------|------------|---------------------|-------------------|
| Snakemake | pull / make (static DAG) | single-lab reproducible research, HPC, Python/pandas-native, file-pattern logic | Workflow Catalog / wrappers | Python shop, cluster data, want the plan visible before an allocation |
| Nextflow | push / reactive dataflow | cloud/production, dynamic pipelines, executor portability | nf-core (largest, curated) | author a scalable/dynamic pipeline, or move one across local/HPC/cloud by profile |
| nf-core | Nextflow, pre-built | running a standard analysis without authoring | the curated catalog itself | a mainstream analysis (rnaseq, sarek, atacseq) already has a community pipeline - adopt, do not build |
| CWL | pull / typed spec | vendor-neutral portability, regulated/clinical, provenance | limited (Dockstore) | hand a pipeline across institutions/platforms, or need a standardized provenance artifact |
| WDL | pull / declarative | GATK/Broad world, NIH-cloud data | WARP (Broad) | controlled-access data in Terra/AnVIL/BioData Catalyst, or GATK Best Practices at scale |

For any mainstream analysis, ADOPT a curated community pipeline (nf-core, WARP) before authoring your own - it already encodes years of QC, edge-case handling, CI tests, and institutional configs. Author from scratch only for genuinely novel logic.

## Sharing and standards

To publish, share, or cite a pipeline so others can run it: Dockstore is the GA4GH Tool Registry Service registry that hosts and versions CWL/WDL/Nextflow/Galaxy workflows (search it before authoring); BioContainers (quay.io) and Bioconda are the standard per-tool container/environment sources engines pull; RO-Crate/CWLProv package a run's provenance portably; and the GA4GH WES/TRS/DRS APIs let a workflow be registered and executed across platforms without changing it. Provenance (the recorded lineage of every output) is distinct from reproducibility (the ability to regenerate it) - CWL ships the strongest standardized provenance artifact out of the box.

## Example Prompts

- "Create a Snakemake workflow for RNA-seq analysis"
- "Set up a Nextflow pipeline with Docker containers"
- "Run nf-core/rnaseq on my samples with an institutional SLURM config"
- "Add a checkpoint for outputs that are unknown until a split step runs"
- "Write a portable CWL workflow and declare the BAM index as a secondaryFile"
- "Create a WDL workflow for GATK variant calling on Terra with dynamic disk sizing"

## Requirements

```bash
# Snakemake (+ HPC executor / cloud storage plugins as needed)
pip install snakemake snakemake-executor-plugin-slurm

# Nextflow (needs a JVM) and the nf-core tools
curl -s https://get.nextflow.io | bash
pip install nf-core

# CWL reference runner
pip install cwltool

# WDL: miniwdl for local dev/linting; Cromwell jar for production/Terra
pip install miniwdl
```

## Related Skills

- **workflows** - End-to-end domain pipelines these engines orchestrate (rnaseq-to-de, fastq-to-variants, and more)
- **read-qc** - QC steps a pipeline wraps as early stages
- **read-alignment** - Alignment steps invoked inside pipeline rules/processes
- **differential-expression** - Downstream analysis steps in an RNA-seq pipeline
