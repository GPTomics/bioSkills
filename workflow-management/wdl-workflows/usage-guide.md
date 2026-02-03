# WDL Workflows - Usage Guide

## Overview

Workflow Description Language (WDL) is a workflow specification designed for data-intensive scientific analysis, widely used for GATK pipelines and cloud execution on Terra/AnVIL platforms.

## Prerequisites
```bash
# Install miniwdl (lightweight runner)
pip install miniwdl

# Or install Cromwell (full-featured)
# Download from: https://github.com/broadinstitute/cromwell/releases
wget https://github.com/broadinstitute/cromwell/releases/download/91/cromwell-91.jar

# Install womtool for validation
wget https://github.com/broadinstitute/cromwell/releases/download/91/womtool-91.jar

# Docker required for containerized execution
```

## Quick Start

Tell your AI agent what you want to do:
- "Create a WDL workflow for variant calling with GATK"
- "Convert this pipeline to WDL for Terra"
- "Run this WDL workflow locally with miniwdl"

## Example Prompts

### Basic Tasks
> "Create a WDL task for BWA-MEM alignment with proper runtime resources"

> "Write a WDL task definition for samtools markdup"

### Workflows
> "Build a WDL workflow for paired-end RNA-seq quantification"

> "Create a scatter workflow to process multiple samples in parallel"

### Cloud Execution
> "Configure my WDL workflow for Google Cloud execution"

> "Set up preemptible instances to reduce costs on Terra"

### GATK Pipelines
> "Create a germline variant calling workflow following GATK best practices"

> "Build a somatic variant calling pipeline with Mutect2"

### Advanced
> "Add struct types to organize my sample inputs"

> "Import subworkflows to modularize my pipeline"

## What the Agent Will Do
1. Create WDL task definitions with inputs, outputs, command, and runtime
2. Wire tasks together in a workflow with proper data flow
3. Add Docker containers for each task
4. Configure runtime resources (CPU, memory, disk)
5. Generate an inputs JSON template
6. Provide commands for validation and execution

## Key Concepts

| Concept | Description |
|---------|-------------|
| Task | Single command/tool definition |
| Workflow | Multiple tasks with data flow |
| Scatter | Parallel execution over arrays |
| Struct | Custom data types for grouped inputs |

## Run Commands
```bash
# Validate WDL syntax
java -jar womtool.jar validate workflow.wdl

# Generate inputs template
java -jar womtool.jar inputs workflow.wdl > inputs.json

# Run with Cromwell
java -jar cromwell.jar run workflow.wdl -i inputs.json

# Run with miniwdl
miniwdl run workflow.wdl -i inputs.json

# Run with miniwdl (specific inputs)
miniwdl run workflow.wdl sample_id=test fastq=reads.fq.gz
```

## Tips
- Always validate WDL with womtool before running
- Use `version 1.0` at the top of every WDL file
- Calculate disk size dynamically based on input file sizes
- Use preemptible instances (preemptible: 3) for cost savings on cloud
- Use structs to organize related inputs (sample ID + files)
- miniwdl is faster for local testing; Cromwell for production/cloud
- GATK workflows are available at github.com/gatk-workflows

## Related Skills
- workflow-management/cwl-workflows - Alternative portable workflow language
- workflow-management/snakemake-workflows - Python-based alternative
- workflow-management/nextflow-pipelines - Groovy-based alternative with nf-core
