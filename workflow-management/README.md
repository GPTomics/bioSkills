# workflow-management

## Overview

Reproducible pipeline frameworks for scalable bioinformatics analyses with dependency management and cluster execution.

**Tool type:** mixed | **Primary tools:** Snakemake, Nextflow, cwltool, Cromwell

## Skills

| Skill | Description |
|-------|-------------|
| snakemake-workflows | Build reproducible pipelines with Snakemake rules and DAGs |
| nextflow-pipelines | Create containerized workflows with Nextflow DSL2 |
| cwl-workflows | Create portable, standards-based pipelines with Common Workflow Language |
| wdl-workflows | Build workflows with WDL for Terra/AnVIL and GATK pipelines |

## Example Prompts

- "Create a Snakemake workflow for RNA-seq analysis"
- "Set up a Nextflow pipeline with Docker containers"
- "Run my workflow on a SLURM cluster"
- "Add checkpointing to my pipeline"
- "Convert my pipeline to CWL for portability"
- "Create a WDL workflow for GATK variant calling"

## Requirements

```bash
# Snakemake
pip install snakemake

# Nextflow
curl -s https://get.nextflow.io | bash

# CWL
pip install cwltool

# WDL
pip install miniwdl
```

## Related Skills

- **workflows** - End-to-end analysis pipelines
- **read-qc** - QC steps in pipelines
- **differential-expression** - Analysis steps
