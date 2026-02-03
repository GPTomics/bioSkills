# CWL Workflows - Usage Guide

## Overview

Common Workflow Language (CWL) is a specification for describing analysis workflows that are portable across execution platforms, making it ideal for sharing pipelines and running on diverse infrastructure.

## Prerequisites
```bash
# Install cwltool (reference implementation)
pip install cwltool

# Or via conda
conda install -c conda-forge cwltool

# For HPC/cloud execution, install Toil
pip install toil[cwl]

# Docker or Singularity for containerized execution
```

## Quick Start

Tell your AI agent what you want to do:
- "Create a CWL workflow for variant calling"
- "Convert my shell pipeline to CWL"
- "Run this CWL workflow with Singularity"

## Example Prompts

### Basic Tools
> "Create a CWL CommandLineTool for running BWA-MEM alignment"

> "Write a CWL tool definition for samtools sort"

### Workflows
> "Build a CWL workflow that runs FastQC, trimming, and alignment"

> "Create a scatter workflow to process multiple samples in parallel"

### Portability
> "Add Docker containers to my CWL tools for reproducibility"

> "Make my workflow compatible with both Docker and Singularity"

### Execution
> "Run this CWL workflow on my local machine"

> "Set up Toil to run my CWL pipeline on a SLURM cluster"

### Advanced
> "Add conditional execution to skip QC if already done"

> "Create a subworkflow for the alignment steps"

## What the Agent Will Do
1. Create CWL tool definitions (CommandLineTool) for each step
2. Wire tools together in a Workflow with proper input/output connections
3. Add container requirements for reproducibility
4. Configure resource requirements (CPU, memory)
5. Create a job input file (YAML) for running the workflow
6. Provide commands for validation and execution

## Key Concepts

| Concept | Description |
|---------|-------------|
| CommandLineTool | Single tool/command definition |
| Workflow | Multiple steps connected by data flow |
| Scatter | Parallel execution over arrays |
| secondaryFiles | Associated files (e.g., .bai with .bam) |

## Run Commands
```bash
# Validate syntax
cwltool --validate workflow.cwl

# Run locally
cwltool workflow.cwl inputs.yaml

# Run with Docker
cwltool --docker workflow.cwl inputs.yaml

# Run with Singularity
cwltool --singularity workflow.cwl inputs.yaml

# Run on HPC with Toil
toil-cwl-runner --batchSystem slurm workflow.cwl inputs.yaml
```

## Tips
- Always validate CWL files before running with `cwltool --validate`
- Use `secondaryFiles` for index files that accompany data files
- The `scatter` feature parallelizes execution over file arrays
- CWL v1.2 supports conditional execution with `when`
- cwltool is the reference implementation; use Toil for HPC/cloud
- Register workflows on Dockstore or WorkflowHub for sharing

## Related Skills
- workflow-management/wdl-workflows - Alternative portable workflow language
- workflow-management/snakemake-workflows - Python-based alternative
- workflow-management/nextflow-pipelines - Groovy-based alternative with nf-core community
