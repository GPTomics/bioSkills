# Nextflow Pipelines Usage Guide

## Overview

Nextflow is a workflow framework for scalable, reproducible pipelines with native container support and cloud execution.

## Quick Start Prompts

- "Create a Nextflow pipeline for RNA-seq analysis"
- "Add Docker containers to my Nextflow processes"
- "Run my pipeline on a SLURM cluster"
- "Set up AWS Batch execution for my workflow"

## Key Concepts

| Concept | Description |
|---------|-------------|
| Process | Execution unit with inputs/outputs |
| Channel | Data flow between processes |
| Module | Reusable process definition |
| Profile | Execution configuration |

## Workflow

1. **Define processes** with containers
2. **Connect** via channels
3. **Configure** execution profile
4. **Run** with `nextflow run`

## Requirements

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Or via conda
conda install -c bioconda nextflow

# Docker or Singularity for containers
```

## Run Commands

```bash
# Basic run
nextflow run main.nf

# With profile
nextflow run main.nf -profile docker

# Resume failed run
nextflow run main.nf -resume

# Run nf-core pipeline
nextflow run nf-core/rnaseq -profile docker
```

## Key Features

- **Container-first** - Docker/Singularity built-in
- **Cloud-native** - AWS, Google, Azure support
- **Scalable** - From laptop to HPC cluster
- **nf-core** - Community pipelines

## Related Skills

- **workflow-management/snakemake-workflows** - Alternative
- **workflows/** - Pre-built analysis pipelines
