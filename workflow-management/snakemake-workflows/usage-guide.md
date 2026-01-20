# Snakemake Workflows Usage Guide

## Overview

Snakemake is a Python-based workflow management system that uses rules to define analysis steps with automatic dependency resolution.

## Quick Start Prompts

- "Create a Snakemake workflow for variant calling"
- "Add cluster execution to my Snakemake pipeline"
- "Set up conda environments for each rule"
- "Generate a DAG visualization of my workflow"

## Key Concepts

| Concept | Description |
|---------|-------------|
| Rule | Single processing step |
| Wildcard | Variable in file patterns |
| Expand | Generate file lists |
| DAG | Directed acyclic graph of dependencies |

## Workflow

1. **Define rules** with input/output patterns
2. **Add target** in `rule all`
3. **Dry run** with `snakemake -n`
4. **Execute** with `snakemake --cores N`

## Requirements

```bash
pip install snakemake

# With conda integration
pip install snakemake pulp

# For cluster execution
pip install snakemake-executor-plugin-slurm
```

## Key Features

- **Automatic parallelization** based on dependencies
- **Conda/container** integration for reproducibility
- **Cluster execution** with profiles
- **Checkpointing** for dynamic workflows
- **Benchmarking** built-in

## Related Skills

- **workflow-management/nextflow-pipelines** - Alternative framework
- **workflows/** - Pre-built analysis pipelines
