# CWL Workflows - Usage Guide

## Overview

Common Workflow Language (CWL) is a SPECIFICATION, not an engine - the fact from which every CWL decision follows. It deliberately splits the portable, strongly-typed workflow DESCRIPTION from its EXECUTION, so one unchanged document runs identically on any conforming runner (cwltool, Toil, Arvados, Calrissian). The entire value proposition - portability, static type-checking, vendor-neutrality, and standardized provenance - is a consequence of that split, and it is why CWL is the most verbose and most explicit of the four major systems: the verbosity buys analyzability and portability. The type system (File/Directory/record/enum/optional, plus secondaryFiles for index companions) catches wiring errors before any compute, and `cwltool --provenance` emits a CWLProv Research Object that hands an auditor exactly what ran - the reason CWL wins in regulated/clinical settings. Adopting CWL buys reproducible workflow LOGIC only: a clean typed DAG over unpinned tools is not reproducible, so containers, references, and seeds must be pinned separately.

## Prerequisites

```bash
# Reference runner (authoring, validation, provenance, local/CI)
pip install cwltool

# Container runtime (one of):
#   Docker 24+   OR   Singularity/Apptainer 3.8+

# At scale, install a production runner instead of cwltool:
pip install "toil[cwl]"          # HPC/cloud batch (Slurm, Kubernetes, AWS)
#   arvados-cwl-runner (Arvados, clinical/enterprise) or Calrissian (Kubernetes) as needed
```

Inputs: one or more `.cwl` documents (CommandLineTool/Workflow) plus a job/input object (`job.yml`) that supplies typed values. Container images pull from BioContainers (quay.io). Target `cwlVersion: v1.2`.

## Quick Start

Tell your AI agent what you want to do:
- "Wrap bwa-mem as a CWL CommandLineTool and declare the reference index as secondaryFiles"
- "Build a CWL workflow that trims with fastp then quantifies with Salmon"
- "Scatter my alignment step over paired FASTQ arrays with dotproduct"
- "Move my ResourceRequirement from hints to requirements so the step stops getting OOM-killed"
- "Validate this workflow and emit a CWLProv provenance object"

## Example Prompts

### Tool Wrapping and Typing
> "Write a CWL CommandLineTool for samtools sort, declaring the input BAM's `.bai` as a secondaryFile and pinning the BioContainers image."

> "My reference FASTA needs a `.fai` and a `.dict` companion - declare both as secondaryFiles and explain the caret rule for the `.dict` name."

### Workflows and Scatter
> "Build a CWL workflow that runs fastp then Salmon, wiring the trimmed reads from the first step into the second with outputSource."

> "Scatter my per-sample alignment over paired R1/R2 arrays - which scatterMethod do I use, and what output shape do I get?"

### Portability and Provenance
> "Audit my workflow for portability leaks: flag every InlineJavascriptRequirement and unpinned container, and rewrite the ${...} expressions as $(...) where possible."

> "Run my workflow with `--provenance` and explain what the CWLProv Research Object contains for an audit."

### Execution
> "Validate this CWL locally with cwltool, then give me the toil-cwl-runner command to run the same document on our SLURM cluster."

## What the Agent Will Do

1. Establish the frame: CWL is a spec; pick the runner (cwltool for authoring/validation, Toil/Arvados/Calrissian for scale) and confirm CWL is the right choice versus Nextflow/WDL/Snakemake.
2. Write CommandLineTool documents that bind typed inputs to the command line and capture outputs via `outputBinding.glob` or `stdout`.
3. Declare `secondaryFiles` on every indexed File input and output (`.bai`, `.fai`/`.dict` with the caret rule, `.tbi`) so runners co-stage companions.
4. Wire tools in a Workflow with explicit `source`/`outputSource` connections and `--validate` the contract before any compute.
5. Place anything that must hold (container, RAM/cores) under `requirements`, not `hints`, and reason about the workflow -> step -> tool override scope.
6. Choose `scatterMethod` deliberately (dotproduct vs flat_/nested_crossproduct) for the intended job cardinality and output shape.
7. Enforce portability discipline: digest-pinned multi-arch containers, `$(...)` over `${...}`, engine extensions kept in hints.
8. Emit a CWLProv provenance object where an audit trail is required, and point to Dockstore/GA4GH TRS for sharing.

## Tips

- CWL is a spec, not cwltool: "CWL is slow" almost always means "I ran cwltool" - move to Toil/Arvados/Calrissian at scale rather than blaming the standard.
- Declare `secondaryFiles` on any indexed File, or the runner stages the primary without its index and the tool fails with "index not found" at runtime.
- The caret `^` strips one extension: `^.dict` on `genome.fasta` yields `genome.dict`, not `genome.fasta.dict`; each leading `^` strips one more.
- Put resources and containers under `requirements` (must hold) not `hints` (advisory) - a `ResourceRequirement` under hints can be silently ignored and OOM-kill the step.
- Requirements inherit Workflow -> step -> tool with the innermost winning; set a default at workflow scope and override per step where a tool needs more.
- Pick `scatterMethod` by intent: `dotproduct` zips equal-length arrays (N jobs); `flat_`/`nested_crossproduct` run all N x M pairs, differing only in output nesting.
- Prefer `$(...)` parameter references (no JS engine, statically analyzable, portable); treat each `${...}` and its `InlineJavascriptRequirement` as a portability debt.
- Pin `DockerRequirement` images by `@sha256:` digest, not `:latest` - a moving tag makes the result non-reproducible no matter how clean the DAG.
- Workflow outputs wire with `outputSource: step/out`; only tool outputs use `outputBinding.glob` - mixing them is a validation error.
- In v1.2, ExpressionTool outputs are not type-checked (a known reference-impl gap); do not lean on ExpressionTool for type safety.
- Use `cwltool --provenance` for a CWLProv Research Object in regulated/audited settings, and register on Dockstore/WorkflowHub (GA4GH TRS) for sharing.

## Related Skills

- workflow-management/wdl-workflows - Alternative portable language for the Terra/GATK ecosystem
- workflow-management/nextflow-pipelines - Reactive-dataflow alternative with the nf-core community catalog
- workflow-management/snakemake-workflows - Python/file-pattern alternative for single-lab HPC
- workflows/fastq-to-variants - An end-to-end variant-calling pipeline that a CWL workflow can orchestrate
