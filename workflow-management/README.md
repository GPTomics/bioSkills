# workflow-management

## Overview

Reproducible pipeline frameworks for scalable bioinformatics analyses with dependency management and cluster execution.

**Tool type:** mixed | **Primary tools:** Snakemake, Nextflow

## Skills

| Skill | Description |
|-------|-------------|
| snakemake-workflows | Build reproducible pipelines with Snakemake rules and DAGs |
| nextflow-pipelines | Create containerized workflows with Nextflow DSL2 |

## Example Prompts

- "Create a Snakemake workflow for RNA-seq analysis"
- "Set up a Nextflow pipeline with Docker containers"
- "Run my workflow on a SLURM cluster"
- "Add checkpointing to my pipeline"

## Requirements

```bash
# Snakemake
pip install snakemake

# Nextflow
curl -s https://get.nextflow.io | bash
```

## Related Skills

- **workflows** - End-to-end analysis pipelines
- **read-qc** - QC steps in pipelines
- **differential-expression** - Analysis steps
