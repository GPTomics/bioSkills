# Snakemake Workflows - Usage Guide

## Overview
Snakemake is a Python-based workflow management system that uses rules to define analysis steps with automatic dependency resolution.

## Prerequisites
```bash
pip install snakemake

# With conda integration
pip install snakemake pulp

# For cluster execution
pip install snakemake-executor-plugin-slurm
```

## Quick Start
Tell your AI agent what you want to do:
- "Create a Snakemake workflow for variant calling"
- "Add cluster execution to my pipeline"
- "Set up conda environments for each rule"

## Example Prompts

### Basic Workflows
> "Create a Snakemake workflow for FASTQ to BAM alignment with BWA"

> "Build a variant calling pipeline using GATK best practices"

### Wildcards and Expansion
> "Set up wildcards to process all samples in my data directory"

> "Use expand() to generate output files for all chromosome/sample combinations"

### Environment Management
> "Add conda environment files for each rule in my workflow"

> "Set up Singularity containers for reproducibility"

### Cluster Execution
> "Configure my Snakemake workflow to run on SLURM cluster"

> "Add resource requirements (memory, threads) to each rule"

### Debugging and Visualization
> "Generate a DAG visualization of my workflow dependencies"

> "Help me debug why this rule isn't running"

## What the Agent Will Do
1. Create a Snakefile with rules for each processing step
2. Define input/output patterns with wildcards for sample handling
3. Add a `rule all` target specifying final outputs
4. Configure resources (threads, memory) for each rule
5. Set up conda/container environments if requested
6. Provide commands for dry-run, execution, and debugging

## Key Concepts

| Concept | Description |
|---------|-------------|
| Rule | Single processing step |
| Wildcard | Variable in file patterns |
| Expand | Generate file lists |
| DAG | Directed acyclic graph of dependencies |

## Tips
- Always do a dry run first with `snakemake -n` to verify the DAG
- Use `--cores all` to automatically use all available cores
- Set `--use-conda` to activate rule-specific environments
- Use `--forcerun rule_name` to re-execute a specific rule
- Add `benchmark: 'benchmarks/{sample}.txt'` to track resource usage
- Use `--rerun-incomplete` to resume after failed runs

## Related Skills
- workflow-management/nextflow-pipelines - Alternative workflow framework with container-first design
- workflows - Pre-built analysis pipelines using Snakemake
