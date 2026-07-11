# Reproducibility by Design - Usage Guide

## Overview

Reproducibility by design treats reproducibility as a property built into a workflow from the first module, not a step added at submission. Its central claim (Cohen-Boulakia 2017) is that availability (code on GitHub with a Docker image) is not reproducibility, and reproducibility is the precondition for reuse: nobody can reuse an analysis they cannot re-run and understand. This skill owns the design and comprehension layer, modular structure with explicit input/output contracts and unambiguous tool identity, and routes the other pillars (environment pinning, testing, provenance) to the skills that cover them.

## Prerequisites

```bash
# Workflow structure reconstruction and visualization (Nextflow)
pip install bioflow-insight

# The engines whose structure this skill reasons about
curl -s https://get.nextflow.io | bash        # Nextflow (needs a JVM)
pip install snakemake                          # Snakemake (--dag / --rulegraph)

# Graphviz to render code-derived DAGs
conda install -c conda-forge graphviz
```

Inputs: an existing pipeline (a `main.nf` or a `Snakefile`), the tool names/versions it invokes, and the paper or README that describes it.

## Quick Start

Tell your AI agent what you want to do:
- "Restructure this monolithic script into modules with explicit input/output contracts"
- "Reconstruct and visualize my Nextflow workflow so I can see the steps"
- "Resolve which tool and version my pipeline actually uses and record it unambiguously"
- "Explain why 'it's on GitHub with Docker' is not the same as reproducible"
- "Decide whether to adopt a curated pipeline instead of maintaining my own"

## Example Prompts

### Designing for Modularity and Reuse
> "Refactor my 300-line analysis script into one-tool-per-step modules with declared inputs and outputs and no hard-coded paths."

> "Review my pipeline for reuse: can each block be understood and run without reading its internals?"

### Comprehension and Structure
> "Run BioFlow-Insight on my main.nf and give me a structure diagram I can keep synced with the code."

> "Generate a code-derived rule graph for my Snakemake workflow instead of my hand-drawn slide."

### Tool Identity and Naming
> "Map the tool names in my paper to the tools my code actually calls, and pin each to a bio.tools ID and container digest."

> "My paper says bwa but the code runs bwa-mem2 at an unstated version, fix the identity record."

### Decision and Routing
> "Is my project reproducible or just available? Tell me which of the five pillars I am missing."

> "Should I keep maintaining my RNA-seq pipeline or adopt nf-core/rnaseq?"

## What the Agent Will Do

1. Name which of the three levels of reproducibility (L1/L2/L3) the request is really about, so the goal is explicit.
2. Separate availability from reproducibility from reuse, and identify which of the five pillars is missing.
3. Restructure logic into modules with explicit input/output contracts, no hard-coded paths, and data separated from code.
4. Establish unambiguous tool identity (bio.tools/EDAM, version, container digest) and reconcile paper-vs-code tool names.
5. Reconstruct the workflow structure from the code (BioFlow-Insight, `--rulegraph`, `-with-dag`) so the diagram cannot drift.
6. Route environment pinning to workflow-management, testing to workflow-testing, and provenance to workflow-provenance.
7. Recommend adopting a curated community pipeline when the analysis is mainstream, since reuse of bespoke modules is empirically rare.

## Tips

- Availability is necessary but is where most projects stop; reproducibility additionally needs the pinned environment, reference data, and enough structure to be understood.
- The unit of reuse is the step, not the pipeline; design each module so a consumer can answer what it does, how to use it, and what it depends on without reading the internals.
- Hard-coded paths mid-script are a top reproducibility killer; parameterize inputs and move data into config/samplesheets.
- A hand-drawn diagram drifts from the code and is worse than none; generate it from the code so it stays true to what runs.
- "Same tool" is ambiguous without a controlled identity; pin a bio.tools ID plus a container digest, not a free-text name.
- Before comparing a new method to an old one, re-run the old pipeline on the old data to confirm the baseline still reproduces.
- For any mainstream analysis, adopt a curated pipeline first; empirically, almost nothing bespoke is reused as-is, while curated workflows are reused much more.

## Related Skills

- reproducibility/workflow-testing - Reference data and regression tests that catch silent breakage
- reproducibility/workflow-provenance - Fine-grained provenance and reproducible publishing
- workflow-management/nextflow-pipelines - Pin containers by digest and use DSL2 module conventions
- workflow-management/snakemake-workflows - Rules, conda lockfiles, and code-derived rule graphs
- workflow-management/nf-core-pipelines - Adopt a curated, reviewed community pipeline
