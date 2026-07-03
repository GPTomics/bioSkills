---
name: bio-workflow-management-nextflow-pipelines
description: Authors reproducible Nextflow DSL2 pipelines built on reactive dataflow, where processes communicate only through channels and execution order is not guaranteed. Use when deciding channel/dataflow (Nextflow) vs rule-based (Snakemake) authoring; wiring queue vs value channels and fixing shared-reference exhaustion with .first(); composing DSL2 modules and subworkflows with take/main/emit; selecting container/conda profiles and pinning images by digest for portability across local/SLURM/LSF/AWS Batch/Google Batch/Kubernetes executors; diagnosing why -resume misses the cache (nondeterministic input order, mtime on network filesystems, mutable :latest tags) with cache 'lenient' and -dump-hashes; managing work/ vs publishDir and dynamic retry escalation; and choosing whether to adopt an nf-core community pipeline or author from scratch.
tool_type: cli
primary_tool: Nextflow
---

## Version Compatibility

Reference examples tested with: Nextflow 24.04+, fastp 0.23+, Salmon 1.10+, MultiQC 1.21+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: Nextflow is DSL2-only (DSL1 was removed in 22.12) and calendar-versioned (24.x/25.x), so any single-script DSL1 tutorial is dead. The `nf-validation` plugin is deprecated in favor of `nf-schema`. A strict, statically-analyzable syntax (VS Code language server) is opt-in now and default in a later release; writing to it future-proofs a pipeline. Pin container images by immutable digest (`@sha256:`), never a moving tag such as `:latest`, or both reproducibility and `-resume` break.

# Nextflow Pipelines

**"Build a scalable, reproducible pipeline with Nextflow"** -> Wire containerized processes together with asynchronous channels (reactive dataflow), so the engine fires each task as soon as its inputs are ready, caches completed tasks for `-resume`, and moves unchanged across executors by swapping a profile.
- CLI: `nextflow run main.nf -profile docker -resume`
- Groovy: DSL2 `process` / `workflow` / channel-operator syntax

## The governing principle: Nextflow is reactive dataflow - processes talk ONLY through channels, execution order is NOT guaranteed

Nextflow is push-model dataflow: processes are pure functions wired together by asynchronous channels, every channel item is a future, and a process fires a task the instant a complete set of inputs is available on ALL its input channels. There is no target file, no backward DAG, no filename-to-rule matching. This is the opposite of Snakemake/CWL/WDL, which are pull/goal-oriented (name a target output, the engine walks a dependency DAG backward to decide what runs). Almost every downstream trap traces back to this one axis:

- The model ENABLES truly dynamic pipelines: a process emits N files computed at runtime, the next process fans out over all N, with no DAG known in advance. Pull engines need a checkpoint/scatter escape hatch for this.
- The COST is that "what will run" is not knowable by reading the script top-to-bottom. Independent branches interleave nondeterministically, so `SAMPLE_A` may finish after `SAMPLE_Z`. Never write logic that assumes order; carry an explicit sample key through the tuple instead.
- "Why did only the first sample run?" and "why did `-resume` re-run everything?" are the two most common support questions, and both are direct symptoms of the dataflow model (see queue-vs-value channels and cache-miss diagnosis below).

A second principle sits above the engine: the DAG buys reproducibility of workflow LOGIC and nothing else automatically (Wratten et al. 2021 *Nat Methods* 18:1161-1168). A clean pipeline over unpinned tools is NOT reproducible. Pin the software environment (container by digest, conda by lockfile), the reference data and params, and control thread/locale/arch leaks (Grüning et al. 2018 *Cell Syst* 6:631-635). The engine gives one layer; the author pins the rest.

## Decision: Nextflow vs the other engines

| Axis | Nextflow (DSL2) | Snakemake | WDL (Cromwell/miniwdl) | CWL |
|------|-----------------|-----------|------------------------|-----|
| Model | reactive dataflow (push, no target) | pull/goal (target -> backward DAG) | pull/goal (declared outputs) | pull/goal (typed, declared) |
| Dynamic DAG (shape depends on runtime data) | native, trivial | `checkpoints` (bolted on) | `scatter` (static-ish) | limited |
| Cloud/executor portability | best-in-class (swap by profile) | good (v8 plugins, catching up) | strong on GCP/Terra | via Toil/Arvados |
| Community pipelines | nf-core (largest, curated) | Workflow Catalog (smaller) | WARP (Broad) | limited |
| Best when | cloud/production, dynamic pipelines, want nf-core | Python shop, HPC, file-pattern logic | Terra/AnVIL, GATK best practices | vendor-neutral portability, regulated |

Honest take: Nextflow wins on executor portability, nf-core, and dynamic pipelines; Snakemake wins on approachability for Python users. Pick by the ecosystem to integrate with, not by benchmarks.

## Decision: queue vs value channel (THE #1 footgun)

| Need | Channel type | Create with | Exhaustion behavior |
|------|--------------|-------------|---------------------|
| One item consumed by one task (per-sample reads) | queue | `Channel.of`, `.fromPath`, `.fromFilePairs`, `.splitCsv` | consumed once, then empty forever |
| A shared value reused on EVERY task (a reference/index) | value (singleton) | `Channel.value(x)`, `.first()`, `.collect()`, or a bare param | read unlimited times, never exhausted |

The firing rule to memorize: a process launches a new task only when EVERY input channel can supply an item; when a queue input drains, no more tasks fire even if other inputs still have items. So a shared reference passed as a queue channel is consumed by the first sample and every later sample silently never runs (exit 0, no error). Corollary: if all of a process's inputs are value channels its outputs are value channels too; if any input is a queue channel the outputs are queue channels.

## Decision: operator selection (and the silent-data-loss traps)

| Operator | Does | Trap / when-wrong |
|----------|------|-------------------|
| `map` | transform each item | pure only; no I/O side effects |
| `collect` | ALL items -> one list item (queue -> value) | blocks until upstream closes; gathers inputs for one aggregating task (MultiQC) |
| `groupTuple` | group by key into `[key, [items]]` | bare form WAITS for the whole channel to close (serialization/deadlock); pass `size: N` or `groupKey(key, n)`; SORT the grouped list or resume breaks |
| `join` | inner-join two channels by key | SILENTLY DROPS non-matching keys by default; use `remainder: true` or `failOnMismatch: true` |
| `combine` | Cartesian product (optional `by:`) | intentional all-vs-all; distinct from `join` (1:1 merge) |
| `mix` | interleave channels into one | order not preserved; pool outputs before a `collect` |
| `branch` | route items to named sub-channels | the DSL2 idiom for conditional routing (single_end vs paired) |
| `first` | first item as a VALUE channel | THE queue -> value fix for shared references |
| `ifEmpty` | supply a default if empty | guards the "empty branch silently vanishes" trap |

## Decision: executor (the portability payoff)

| Target | `executor` | Best when | Watch |
|--------|-----------|-----------|-------|
| Laptop/dev | `local` | development, tiny data, `-stub` wiring tests | one machine only |
| On-prem HPC | `slurm`, `lsf`, `sge`, `pbs` | shared cluster, on-prem data | tune `queueSize`/`submitRateLimit`; `scratch true` on slow shared FS |
| Cloud batch | `awsbatch`, `google-batch`, `azurebatch` | elastic scale, no on-prem HPC | input localization copy dominates cost/time; Wave + Fusion cut it |
| Kubernetes | `k8s` | already running K8s | more setup overhead |

Never bake the executor into pipeline code; always set it in a profile so the same code moves across all of them.

## Process and DSL2 modules (take / main / emit)

A DSL2 module wraps one tool as a `process` and can be `include`d and called multiple times (aliased), which DSL1 could not. Subworkflows compose modules with named inputs/outputs.

```groovy
// modules/fastqc.nf -- one tool, reusable, tested in isolation
process FASTQC {
    tag "${meta.id}"                                   // meta map threads sample identity through every operator
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'  // pin an immutable tag/digest, never :latest
    label 'process_low'                                // maps to central resource config, decoupled from module code

    input:
    tuple val(meta), path(reads)                       // nf-core convention: [ meta, files ], meta = [id:'x', single_end:false]

    output:
    tuple val(meta), path('*.zip'), emit: zip          // meta round-trips so downstream always knows the sample

    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}
```

```groovy
// subworkflows/qc.nf -- take/main/emit names the interface
include { FASTQC } from '../modules/fastqc'
include { MULTIQC } from '../modules/multiqc'

workflow QC {
    take:
    reads

    main:
    FASTQC(reads)
    MULTIQC(FASTQC.out.zip.collect())                  // collect() gathers all samples' zips into ONE aggregating task

    emit:
    report = MULTIQC.out.report                        // access as QC.out.report from the caller
}
```

## Channels: the shared-reference fix, explained

```groovy
workflow {
    reads_ch = Channel.fromFilePairs(params.reads)     // queue: [id, [r1, r2]] per sample -- consumed once each
    index_ch = Channel.fromPath(params.index)          // queue: ONE item, the shared index

    // BUG if written ALIGN(reads_ch, index_ch): the index is consumed by sample 1,
    // its queue is then empty, and samples 2..N silently never fire (exit 0, no error).
    // .first() converts the queue to a VALUE channel, reusable on every task invocation.
    ALIGN(reads_ch, index_ch.first())
}
```

## Resume: the task hash, and diagnosing a cache miss

`-resume` reuses a task only on an EXACT hit of the task hash, computed from the input file identities, the resolved `script` text, the container reference, and input values/params. A single-bit change in any component busts the cache and re-runs the task. The notorious silent causes:

- Nondeterministic input ordering (`collect`/`groupTuple`/glob expansion) -> the ordered list is part of the hash. Fix: `toSortedList()` or `.map{ k, v -> [k, v.sort()] }`.
- Mutable container tags (`:latest`, or a re-pushed version) -> pin by digest.
- mtime hashing on network filesystems (Lustre/NFS) -> `-resume` works locally but misses on the cluster. Fix: `cache 'lenient'` (hashes size + path, ignores mtime).
- Absolute paths, dates, `$RANDOM`, or `hostname` baked into the `script` string -> the script hash changes every run.
- A deleted or moved `work/` -> the hash hits the DB but the task dir is gone, forcing a re-run.

```groovy
process ALIGN {
    // 'lenient' skips mtime -- the single most useful resume fix on HPC/cloud shared filesystems.
    // 'deep' hashes full file CONTENT (slower, robust when metadata lies); 'false' never caches.
    cache 'lenient'
    // ...
}
```

Definitive diagnosis: run both executions with `-dump-hashes` and diff which hash component differed, or `nextflow log <run_name> -f hash,name,status,workdir` to compare per-task hashes across runs. Everything else is guessing.

## work/ is truth; publishDir is a side effect

Every task runs in an isolated `work/<hash>/` dir holding the real outputs plus the forensic trail (`.command.sh` resolved script, `.command.log`, `.exitcode`). That directory IS the pipeline's output store and the ONLY thing `-resume` reads. `publishDir` merely copies or symlinks SELECTED outputs to a human-friendly location, and its failure can be SILENT because the task itself exited 0 in `work/`. Consequences:

- `mode: 'symlink'` (default) breaks if `work/` is later deleted; `mode: 'copy'` is safe to delete afterward; `mode: 'move'` breaks `-resume` (the output leaves `work/`), so use it only for terminal outputs.
- Never `rm -rf work/` if a resume might be wanted; use `nextflow clean` (which prunes the cache DB consistently). "Outputs missing but the pipeline succeeded" almost always means looking in `publishDir` instead of `work/<hash>/`.

## Resources: dynamic retry escalation

```groovy
process BIG {
    // 137=SIGKILL/OOM, 143=SIGTERM (SLURM wall-time kill); the 130..145 signal band + 104 (transient I/O) retry, fail fast otherwise.
    errorStrategy { task.exitStatus in ((130..145) + 104) ? 'retry' : 'terminate' }
    maxRetries 3
    memory { 8.GB * task.attempt }                     // task.attempt is 1-based; escalates 8 -> 16 -> 24 -> 32 GB
    time   { 4.h  * task.attempt }                     // a transient OOM auto-escalates instead of killing the run

    script:
    """
    memory_intensive_command
    """
}
```

`errorStrategy` values are `'terminate'` (default), `'retry'`, `'ignore'` (drop the failed task's outputs and continue over survivors), and `'finish'` (graceful drain). The nf-core `process.resourceLimits` directive (Nextflow 24.04+, which replaced the pre-3.0 `check_max` pattern) clamps the escalated request to the machine/queue ceiling so `8.GB * task.attempt` never asks for more than a node has.

## Executors and profiles: portability

```groovy
// nextflow.config -- executor lives in a profile, never in the pipeline code
profiles {
    docker      { docker.enabled = true }
    singularity { singularity.enabled = true }
    slurm {
        process.executor = 'slurm'
        executor { queueSize = 100; submitRateLimit = '10/1min' }   // avoid hammering the scheduler
    }
    awsbatch {
        process.executor = 'awsbatch'
        aws.region = 'us-east-1'
    }
}

process {
    cpus = 2; memory = '4 GB'; time = '1h'             // sane defaults
    withLabel: 'process_high' { cpus = 16; memory = '64 GB'; time = '12h' }  // labels centralize per-tier tuning
}
```

Run with `-profile slurm,singularity` (comma-separated, NO spaces; later profiles override earlier).

## Adopt an nf-core pipeline before authoring

For any mainstream analysis (RNA-seq, variant calling, ATAC, methylation, amplicon), a curated `nf-core/<pipeline>` already encodes years of QC, containerized modules, nf-test regression tests, and institutional configs. Reinventing it is months of work and worse QC. Pin the revision: `nextflow run nf-core/rnaseq -r 3.14.0 -profile test,docker --outdir results`. DIY is justified only for genuinely novel logic. See workflow-management/nf-core-pipelines for running, configuring, and building samplesheets against community pipelines; this skill covers AUTHORING.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Only the first sample processed, exit 0, no error | shared reference on a queue channel, exhausted after task 1 | `.first()` / `Channel.value` on the reference |
| Pipeline hangs at a grouping step | `groupTuple` with no size on a channel that never closes | `size: N` or `groupKey(key, n)` |
| Some samples silently disappear mid-pipeline | `join` dropped non-matching keys | `remainder: true` or `failOnMismatch: true` |
| `-resume` re-runs everything | nondeterministic input order, or `:latest` tag, or mtime on network FS | sort inputs; pin container digest; `cache 'lenient'` |
| Resume works locally, misses on the cluster | mtime unreliable on Lustre/NFS | `cache 'lenient'` |
| Outputs missing but the pipeline "succeeded" | publishDir failed silently, or looked in publishDir not work/ | check `work/<hash>/`; use `mode: 'copy'` |
| Resume broken after cleanup | deleted `work/` | never `rm -rf work/`; use `nextflow clean` |
| OOM kills a long run near the end | fixed memory, no escalation | `memory { 8.GB * task.attempt }` + conditional retry |
| Wrong result, no error, after a base image update | mutable tag served a stale cache hit | pin by digest; `cache 'deep'` for critical inputs |
| Huge cloud bill / slow S3 pipeline | explicit stage-in/out copies of large files | Wave + Fusion (POSIX over object store) |
| `-profile test docker` ignores docker | space instead of comma | `-profile test,docker` |

## Related Skills

- workflow-management/nf-core-pipelines - Run and configure community pipelines (the RUNNING counterpart to this authoring skill)
- workflow-management/snakemake-workflows - Pull/goal-oriented alternative for Python shops and file-pattern logic
- workflows/rnaseq-to-de - End-to-end RNA-seq quantification to differential expression
- read-qc/quality-reports - QC steps a pipeline orchestrates (FastQC/MultiQC)

## References

- Di Tommaso P, Chatzou M, Floden EW, Prieto Barja P, Palumbo E, Notredame C. 2017. Nextflow enables reproducible computational workflows. *Nat Biotechnol* 35(4):316-319.
- Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. 2020. The nf-core framework for community-curated bioinformatics pipelines. *Nat Biotechnol* 38(3):276-278.
- Wratten L, Wilm A, Göke J. 2021. Reproducible, scalable, and shareable analysis pipelines with bioinformatics workflow managers. *Nat Methods* 18:1161-1168.
- Grüning B, Chilton J, Köster J, et al. 2018. Practical computational reproducibility in the life sciences. *Cell Syst* 6(6):631-635.
