# Nextflow Pipelines - Usage Guide

## Overview
Nextflow is a workflow framework for scalable, reproducible pipelines with native container support and cloud execution.

## Prerequisites
```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Or via conda
conda install -c bioconda nextflow

# Docker or Singularity for containers
```

## Quick Start
Tell your AI agent what you want to do:
- "Create a Nextflow pipeline for RNA-seq analysis"
- "Add Docker containers to my processes"
- "Run my pipeline on a SLURM cluster"

## Example Prompts

### Basic Pipelines
> "Create a Nextflow pipeline for paired-end read alignment with STAR"

> "Build a variant calling workflow with GATK HaplotypeCaller"

### Containers and Reproducibility
> "Add Docker containers to each process in my pipeline"

> "Configure Singularity for HPC execution"

### Channel Operations
> "Set up channels to pair FASTQ files with sample metadata"

> "Use channel operators to split, combine, and filter data"

### Cluster and Cloud Execution
> "Configure my pipeline to run on SLURM with resource management"

> "Set up AWS Batch execution for my workflow"

### Using nf-core
> "Run the nf-core/rnaseq pipeline on my data"

> "Customize nf-core pipeline parameters for my analysis"

## What the Agent Will Do
1. Create a main.nf workflow file with DSL2 syntax
2. Define processes with input/output channels and containers
3. Set up a nextflow.config with execution profiles
4. Configure resource allocation per process
5. Provide commands for local and cluster execution
6. Set up resume capability for failed runs

## Key Concepts

| Concept | Description |
|---------|-------------|
| Process | Execution unit with inputs/outputs |
| Channel | Data flow between processes |
| Module | Reusable process definition |
| Profile | Execution configuration |

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

## Tips
- Always use `-resume` to continue from cached results after failures
- Set `publishDir` in processes to copy outputs to a results directory
- Use `-profile docker,test` to run nf-core pipelines with test data first
- Check execution with `nextflow log` and inspect work directories
- Use `params.outdir` pattern for configurable output locations
- Consider nf-core pipelines before building from scratch

## Related Skills
- workflow-management/snakemake-workflows - Python-based alternative with make-like syntax
- workflows - Pre-built analysis pipelines
