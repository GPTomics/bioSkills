# WDL Workflows - Usage Guide

## Overview

WDL (Workflow Description Language) is chosen for its ecosystem, not its syntax: it is the language of the GATK/Broad and Terra/AnVIL/BioData Catalyst cloud world, hosted on Dockstore and hardened in the WARP pipeline corpus. The governing idea an author must internalize is that on a cloud backend the engine spins up a fresh VM per task, so the `runtime` block is a cost-and-reliability contract - it sizes CPU, memory, and disk, buys preemptible VMs, and sets retries. Localization (copying every input `File` from object storage onto the VM before the command runs) is the dominant cost and wall-time factor, which is why dynamic disk sizing, call caching, and preemptibles all exist. A clean WDL over an unpinned `:latest` container is not reproducible: the engine pins step order, but the author must pin the container by `@sha256:` digest, the reference build, and the parameters.

## Prerequisites

```bash
# miniwdl: lightweight Python engine + static linter (local dev, CI, debugging)
pip install miniwdl

# Cromwell + womtool: the JVM production engine (Terra backend) and its validation toolkit
# Download the matching release jars from github.com/broadinstitute/cromwell/releases
#   cromwell-<ver>.jar, womtool-<ver>.jar   (both require a JVM: java -version)

# Docker (or Podman/Singularity for miniwdl) is required for containerized task execution.
```

Inputs: a `.wdl` workflow file, an `inputs.json` with fully namespaced keys, and container images (pinned by digest). Curated production pipelines to imitate or run: WARP (github.com/broadinstitute/warp) and Dockstore (dockstore.org).

## Quick Start

Tell your AI agent what you want to do:
- "Create a WDL workflow for GATK germline variant calling for Terra"
- "Add dynamic disk sizing to my WDL alignment task so it stops running out of space"
- "Set up preemptible VMs on the short scatter tasks to cut cloud cost"
- "Lint this WDL with miniwdl and tell me why Cromwell won't run it"
- "Convert my parallel-array scatter to a struct so samples stop getting mismatched"

## Example Prompts

### Ecosystem targeting
> "I have controlled-access data in AnVIL. Write a WDL workflow I can run on Terra, and tell me which parts Terra manages for me versus what I configure."

### Runtime as a cost contract
> "Add a `runtime` block to each task that sizes disk dynamically from input file size, pins Docker by digest, and uses preemptible VMs only on the tasks short enough to justify them."

### Scatter and structs
> "Refactor my workflow that scatters over three parallel arrays (ids, R1, R2) into a scatter over an Array of structs so a sample can't get a mismatched read pair."

### Call caching
> "My identical Cromwell reruns recompute everything. Explain what call caching needs to actually cache, and which of my settings are silently busting it."

### Debugging
> "Cromwell gave me a NullPointerException stack trace. Show me how to reproduce the failure under miniwdl to get a readable error."

## What the Agent Will Do

1. Confirm the target ecosystem (Terra/AnVIL/GATK) and whether a WARP/Dockstore pipeline already exists before authoring from scratch.
2. Declare each `task` with typed inputs, a `command <<< >>>` heredoc using `~{}` placeholders, typed outputs, and a `runtime` block.
3. Size runtime disk dynamically from `size()` with a multiplier plus headroom, and pin every `docker:` by `@sha256:` digest.
4. Wire tasks in a `workflow`, using `scatter` over an `Array` (or `Array` of structs) for parallel fan-out.
5. Choose preemptible vs on-demand per task by length and idempotency, and set `maxRetries`/`bootDiskSizeGb` where needed.
6. Lint with `miniwdl check`, validate and scaffold inputs with `womtool`, then run locally (miniwdl) or at scale (Cromwell/Terra).

## Structs, subworkflows, and input JSON namespacing

Bundle related per-sample fields into a `struct` and scatter over an array of them so a read pair cannot be mismatched by parallel-array indexing:

```wdl
version 1.0

struct SampleFastqs { String sample_id; File fastq_1; File fastq_2 }

workflow paired_alignment {
    input { Array[SampleFastqs] samples; File reference }
    scatter (s in samples) {
        call align { input: sample_id = s.sample_id, reads_1 = s.fastq_1, reads_2 = s.fastq_2, reference = reference }
    }
    output { Array[File] bams = align.bam }
}
```

Compose larger pipelines from `import`ed subworkflows, aliasing to avoid namespace clashes:

```wdl
import "qc.wdl" as qc
import "align.wdl" as align

call qc.quality_control { input: reads_1 = fastq_1, reads_2 = fastq_2 }
call align.alignment { input: reads_1 = quality_control.trimmed_1, reference = reference }
```

Input JSON keys are the fully namespaced dotted path `Workflow.[subworkflow.]call_alias.input_name`; optional inputs may be omitted. Generate the skeleton with `womtool inputs`:

```json
{
    "paired_alignment.reference": "gs://bucket/genome.fa",
    "paired_alignment.samples": [
        {"sample_id": "s1", "fastq_1": "gs://bucket/s1_R1.fq.gz", "fastq_2": "gs://bucket/s1_R2.fq.gz"}
    ]
}
```

## Engine CLI cheat

```bash
miniwdl check workflow.wdl              # static lint + ShellCheck; --strict gates CI; # !WarningName suppresses
womtool validate workflow.wdl          # Cromwell-side structural validation
womtool inputs workflow.wdl > inputs.json   # scaffold the namespaced input JSON
womtool graph workflow.wdl             # emit the DAG as DOT
miniwdl run workflow.wdl -i inputs.json     # local run, readable errors
miniwdl run workflow.wdl sample_id=test reads_1=r1.fq.gz   # inputs on the command line
java -jar cromwell.jar run workflow.wdl -i inputs.json      # one-off local/HPC run
java -jar cromwell.jar server           # REST server - the mode Terra uses
```

## Tips

- Size disk dynamically: `ceil(size(f,"GiB")*factor)+buffer`. A static `disks: "local-disk 100 HDD"` fails late (after localization) on the one sample bigger than guessed.
- Use `size(f, "GiB")` consistently (binary GiB); `"GB"` is decimal. Always `ceil()`, never `round()` - disk must round UP.
- Pin every container by `@sha256:` digest, not `:latest` or a mutable tag: this fixes reproducibility AND stabilizes call caching in one move.
- `preemptible: 3` is an Int (number of preemptible attempts before on-demand fallback), not a Boolean. Use preemptibles for short (<~2-4h), idempotent, restart-safe tasks; on-demand for long or stateful ones.
- Use `~{}` inside `command <<< >>>`, never `${}` - the latter collides with bash parameter expansion and is a classic silent bug.
- Call caching is NOT on by default: it needs `call-caching.enabled=true` plus a persistent MySQL/Postgres DB (the default HSQLDB loses the cache on restart). Terra manages this for you.
- When Cromwell throws a JVM stack trace, reproduce the run under `miniwdl run` for an error that points at the WDL line and the failing shell command.
- Do not reinvent joint genotyping or CRAM->GVCF: search WARP/Dockstore/GATK first; those pipelines already encode the disk-sizing and preemptible lessons.
- Remember `Directory` is a 1.2 feature (not 1.1), and `min`/`max`/`None` are 1.1: set the `version` header to what your engine supports.

## Related Skills

- workflow-management/cwl-workflows - Vendor-neutral portable spec; prefer it for cross-institution portability
- workflow-management/nextflow-pipelines - Channel/dataflow engine + nf-core for dynamic/streaming cloud pipelines
- workflow-management/snakemake-workflows - Python-native pull engine for HPC and file-pattern logic
- workflows/fastq-to-variants - The end-to-end variant-calling analysis a GATK WDL orchestrates
- variant-calling/gatk-variant-calling - The GATK Best Practices steps WDL encodes for Terra/WARP
