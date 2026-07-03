---
name: bio-workflow-management-wdl-workflows
description: Authors bioinformatics pipelines in WDL (Workflow Description Language) run by Cromwell or miniwdl, targeting the GATK/Broad and Terra/AnVIL/BioData Catalyst cloud ecosystem, with tasks, workflows, scatter-gather parallelism, structs, and a runtime block that sizes the cloud VM. Use when deciding to target Terra/AnVIL/GATK/WARP (chosen for the ecosystem, not the language); sizing runtime disks dynamically for a fresh-per-task cloud VM (ceil(size(f)*factor)+buffer); choosing preemptible vs on-demand VMs by task length and idempotency; picking Cromwell (production, cloud, call-caching) vs miniwdl (local dev, miniwdl check linting, readable errors); enabling and debugging call-caching silent-miss modes; pinning Docker by digest for reproducibility and cache stability; or scattering an array for parallel fan-out.
tool_type: cli
primary_tool: cromwell
---

## Version Compatibility

Reference examples tested with: Cromwell 87+, miniwdl 1.12+, WDL spec 1.0/1.1/1.2

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: every WDL file must open with a `version 1.0`/`1.1`/`1.2` header; omitting it selects the old draft-2 dialect with no `~{}` interpolation. The first-class `Directory` type is a 1.2 feature, NOT 1.1; `min`/`max`/`None` arrived in 1.1. Cromwell is the JVM production engine that powers Terra; miniwdl is the Python engine used for local dev, static linting (`miniwdl check`), and readable errors. Pin every `docker:` by `@sha256:` digest, never a floating tag.

# WDL Workflows

**"Build a WDL pipeline for Terra/AnVIL or a GATK best-practices run"** -> Declare `task`s (a containerized command with typed inputs/outputs and a runtime block) and wire them in a `workflow`, then run on Cromwell (cloud/Terra) or miniwdl (local).
- CLI: `womtool validate` / `womtool inputs` (Cromwell toolkit), `miniwdl check` (static lint + ShellCheck), `cromwell run` / `miniwdl run` (execute)
- WDL: `version` header, `task`/`workflow`/`call`, `scatter` fan-out, `runtime { docker, cpu, memory, disks }`

## The governing principle: the runtime block is a cost + reliability CONTRACT

WDL is not chosen on language merits; it is the language of a gravitational system - GATK Best Practices -> Cromwell -> Terra/AnVIL/BioData Catalyst -> Dockstore -> WARP (Van der Auwera & O'Connor 2020). One targets WDL because the data or the collaborators already live in that NIH-cloud ecosystem, and to run vetted GATK pipelines without reinventing them. The design bet is human readability over expressive power. The corollary that governs every real decision: on a cloud backend the engine spins up a FRESH VM per task, so the author must declare its CPU, memory, and disk. The `runtime` block is therefore a cost-and-reliability contract, not decoration, and three traps follow from it:

- Localization dominates cost and wall-time. The engine COPIES (localizes) every input `File` from object storage onto the VM's local disk before the command runs, then delocalizes outputs back. A 30 GB CRAM's transfer can dwarf the compute. Disk math, call caching, and preemptibles all exist to manage bytes moved - subset early and avoid re-localizing the same reference into every scatter shard.
- Under-sized `disks` kills the job LATE. A static `disks: "local-disk 100 HDD"` fails on the one sample bigger than guessed, after an hour of localization, with a cryptic "No space left on device". Size disk dynamically from `size()`.
- A pipeline without pinned containers is not reproducible. `docker: "gatk:latest"` silently breaks reproducibility AND busts call caching, because the cache key hashes the resolved image identity. Pin by `@sha256:` digest.

A clean WDL over unpinned tools is not reproducible: the engine pins step order (layer 1); the author must still pin the container by digest, the reference build, and the parameters.

## Decision: choose WDL, and choose its engine

| Author picks WDL when... | Fails / friction when... |
|--------------------------|--------------------------|
| Controlled-access data is in AnVIL/Terra/BioData Catalyst | The pipeline is dynamic/streaming (WDL has no channels; use nextflow-pipelines) |
| Running GATK Best Practices at population scale | Maximum vendor-neutral portability across institutions is the goal (use cwl-workflows) |
| A WARP/Dockstore pipeline already encodes the analysis | Tight Python/pandas HPC integration is wanted (use snakemake-workflows) |

| Engine | Runtime | Reach for it when | Weakness |
|--------|---------|-------------------|----------|
| Cromwell | Scala/JVM | Production cloud, Terra, robust call caching at scale | Cryptic JVM errors; needs MySQL/Postgres for persistent cache; slow startup |
| miniwdl | Python | Local dev, CI, debugging; `miniwdl check` static lint + ShellCheck; readable errors | Not the Terra engine; smaller cloud story |
| womtool | JVM utility | `validate`, generate the inputs JSON skeleton, `graph` the DAG | Not an executor - validation only |

Practical loop: author and lint with `miniwdl check` locally -> validate and scaffold inputs with `womtool` -> run at scale on Cromwell/Terra. When Cromwell throws a JVM stack trace, reproduce under `miniwdl run` for a message that points at the WDL line.

| Factor | Preemptible / spot (`preemptible: N`) | On-demand |
|--------|---------------------------------------|-----------|
| Cost | ~60-91% cheaper | full price |
| Interruption | reclaimable any second, work discarded | stable |
| Fit | short (<~2-4h), idempotent, restart-safe, scatter shards | long, stateful, near-deadline, non-idempotent |
| Anti-pattern | long non-idempotent task -> retry thrash, can cost MORE than on-demand | over-paying for a trivially restartable 20-min task |

`preemptible: 3` is an Int (retry on a preemptible VM up to 3 times, then fall back to on-demand), NOT a Boolean.

## Task and workflow: the reference shape

A `task` bundles a container, typed inputs, a heredoc command with `~{}` placeholders, typed outputs, and a runtime block. A `workflow` `call`s tasks and passes one call's output to the next by name.

```wdl
version 1.0

task fastp {
    input {
        String sample_id
        File reads_1
        File reads_2
        Int threads = 4
    }
    # ~{} is the WDL-idiomatic placeholder; ${} collides with bash parameter expansion.
    command <<<
        fastp -i ~{reads_1} -I ~{reads_2} \
            -o ~{sample_id}_R1.fq.gz -O ~{sample_id}_R2.fq.gz \
            --json ~{sample_id}.json --thread ~{threads}
    >>>
    output {
        File trimmed_1 = "~{sample_id}_R1.fq.gz"
        File trimmed_2 = "~{sample_id}_R2.fq.gz"
    }
    runtime {
        docker: "quay.io/biocontainers/fastp@sha256:<digest>"   # digest, not :latest
        cpu: threads
        memory: "4 GB"
    }
}

workflow trim {
    input { String sample_id; File r1; File r2 }
    call fastp { input: sample_id = sample_id, reads_1 = r1, reads_2 = r2 }
    output { File out_1 = fastp.trimmed_1 }
}
```

## Scatter: explicit parallel fan-out (no channels)

WDL parallelism is explicit: build an `Array`, `scatter` over it (implicitly parallel), and the engine auto-gathers each shard's output into an `Array` in input order. There is no lazy channel to drain.

```wdl
scatter (idx in range(length(sample_ids))) {
    call align {
        input: sample_id = sample_ids[idx], reads = fastq_files[idx], reference = reference
    }
}
# align.bam outside the scatter is an Array[File], gathered in input order.
output { Array[File] bams = align.bam }
```

Bundle per-sample fields into a `struct` (`struct SampleData { String id; File bam }`) and scatter over `Array[SampleData]` to avoid parallel-array index bugs; see usage-guide.md.

## Runtime as a cost contract + dynamic disk sizing

**Goal:** Size the fresh cloud VM so the task neither fails on disk nor over-pays.

**Approach:** Compute disk from actual input size with a multiplier for outputs/intermediates plus headroom, round UP with `ceil()`, and make it overridable.

```wdl
task bwa_mem {
    input { File reads_1; File reads_2; File reference; Int? override_disk_gb }
    # size(f,"GiB") is binary GiB (be consistent); *2.5 covers input+output+intermediates,
    # +20 is headroom. ceil() always rounds UP - disk must never under-size.
    Int disk_gb = select_first([override_disk_gb,
                  ceil((size(reads_1, "GiB") + size(reads_2, "GiB") + size(reference, "GiB")) * 2.5) + 20])
    command <<< bwa mem ~{reference} ~{reads_1} ~{reads_2} > aligned.sam >>>
    output { File sam = "aligned.sam" }
    runtime {
        docker: "quay.io/biocontainers/bwa@sha256:<digest>"
        cpu: 8
        memory: "16 GB"
        disks: "local-disk ~{disk_gb} HDD"   # mount, GB Int, type; HDD cheap/slow, SSD fast/pricey
        bootDiskSizeGb: 20                    # boot disk holds the image; raise for large images
        preemptible: 3                        # Int = # attempts, then on-demand fallback
        maxRetries: 1                         # retries on ANY failure (distinct from preemptible)
    }
}
```

## Call caching: the `-resume` analog, and how it silently misses

Cromwell hashes each call from its command template, input values (including file CONTENT hashes), Docker image identity, and runtime attributes; on a rerun an identical hash reuses prior outputs. Unlike Nextflow's `-resume`, it is NOT on by default: the in-memory HSQLDB loses the cache on restart, so a persistent DB plus config is required (Terra manages this behind a checkbox).

```
call-caching { enabled = true, invalidate-bad-cache-results = true }
# plus a MySQL/PostgreSQL database stanza - the default HSQLDB does not persist the cache.
```

Silent-miss modes: a floating `:latest` tag resolves to a new digest -> new hash -> miss (pin by digest); a touched/re-staged input whose content or mtime changed busts the cache; a path-based hashing strategy misconfigured on a container backend disables caching; and any whitespace change in the `command` block changes the hash.

## Validate, generate inputs, run

```bash
miniwdl check workflow.wdl              # static lint + ShellCheck (add --strict to gate CI)
womtool validate workflow.wdl          # Cromwell-side structural validation
womtool inputs workflow.wdl > inputs.json   # scaffold the namespaced input JSON
miniwdl run workflow.wdl -i inputs.json     # local run, readable errors
java -jar cromwell.jar run workflow.wdl -i inputs.json   # one-off; `cromwell server` = REST (Terra mode)
```

Input JSON keys are fully namespaced `Workflow.[subworkflow.]call_alias.input_name`, e.g. `{"rnaseq.fastp.threads": 8}`; optional inputs may be omitted. See usage-guide.md for structs, subworkflows, and the full namespacing rules.

Do not hand-roll joint genotyping or CRAM->GVCF: WARP publishes production-vetted, cost-tuned WDL to imitate for disk and preemptible discipline (WARP team 2025).

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Task dies late with "No space left on device" | static under-sized `disks` | dynamic `ceil(size(f,"GiB")*factor)+buffer` |
| Job cost balloons; wall-time is mostly "waiting" | localizing huge inputs to every scatter shard | subset early; co-locate data + compute zones; reuse the reference where the backend caches it |
| Reruns recompute everything | call caching off, no persistent DB, or a floating docker tag | enable caching + MySQL/Postgres + digest-pin docker |
| Preemptible task never finishes, costs more than on-demand | long non-idempotent task on `preemptible: N` | move to on-demand or shorten/checkpoint the task |
| VM fails to boot | Docker image larger than the boot disk | raise `bootDiskSizeGb` |
| "Works on my Cromwell, not on Terra" | env drift not baked into the container; unpinned tag | bake everything into a digest-pinned image |
| Cryptic JVM stack trace from Cromwell | engine surfacing an internal error | reproduce under `miniwdl run` for a legible, line-pointing message |
| `${VAR}` in a command expands wrong or breaks | `${}` collides with bash parameter expansion | use `~{}` for WDL interpolation inside `command <<< >>>` |
| Engine rejects `Directory` under `version 1.1` | first-class `Directory` is a 1.2 feature | move the header to `version 1.2` (or `version development` on old engines) |

## Related Skills

- workflow-management/cwl-workflows - Vendor-neutral portable spec; choose it over WDL when handing a pipeline across institutions
- workflow-management/nextflow-pipelines - Channel/dataflow engine + nf-core; choose it for dynamic/streaming cloud pipelines
- workflow-management/snakemake-workflows - Python-native pull engine for HPC and file-pattern logic
- workflows/fastq-to-variants - The end-to-end variant-calling analysis a GATK WDL orchestrates
- variant-calling/gatk-variant-calling - The GATK Best Practices steps WDL encodes for Terra/WARP

## References

- Van der Auwera GA, Carneiro MO, Hartl C, et al. 2013. From FastQ data to high-confidence variant calls: the Genome Analysis Toolkit best practices pipeline. *Curr Protoc Bioinformatics* 43:11.10.1-11.10.33.
- Voss K, Gentry J, Van der Auwera G. 2017. Full-stack genomics pipelining with GATK4 + WDL + Cromwell. *F1000Research* 6:1379 (ISCB Comm J, poster).
- Van der Auwera GA, O'Connor BD. 2020. *Genomics in the Cloud: Using Docker, GATK, and WDL in Terra.* O'Reilly Media. ISBN 9781491975190.
- WARP team (Broad Institute). 2025. WARP analysis research pipelines: cloud-optimized workflows for biological data processing and reproducible analysis. *Bioinformatics* 41(10):btaf494.
- OpenWDL specification. github.com/openwdl/wdl - versioned SPEC on branches wdl-1.0, wdl-1.1 (1.1.3), wdl-1.2; docs at docs.openwdl.org.
