---
name: bio-workflow-management-snakemake-workflows
description: Authors reproducible bioinformatics pipelines with Snakemake - rules wired by output-file pattern, wildcards and expand() for sample fan-out, checkpoints for runtime-unknown outputs, resource/retry escalation, and conda/container software deployment on HPC and cloud. Use when deciding rule-based (Snakemake) vs channel/dataflow (Nextflow) authoring; wiring rules by OUTPUT-file pattern rather than imperative order; using wildcards + expand() for sample fan-out and constraining them to stop silent mis-routing; adding checkpoints when the set of outputs is unknown until a step runs (dynamic DAG); diagnosing why a job reran (or did not) under the mtime-plus-provenance trigger set; escalating memory on retry for OOM-killed jobs; and porting a Snakemake 7 `--cluster`/remote-provider command to the Snakemake 8+ executor-plugin and storage-plugin model (snakemake-executor-plugin-slurm) with `--software-deployment-method`.
tool_type: python
primary_tool: Snakemake
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: Snakemake 8.0+, Python 3.11+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: Snakemake 8 (Jan 2024) removed `--cluster`, `--drmaa`, and the `*RemoteProvider` classes from core and moved them to pip-installable EXECUTOR plugins (`--executor slurm`, package `snakemake-executor-plugin-slurm`) and STORAGE plugins (`storage.s3(...)`, `snakemake-storage-plugin-s3`). `--use-conda`/`--use-singularity` became `--software-deployment-method` / `--sdm conda apptainer`. A Snakemake 7 command line does not run unchanged on 8/9. Run `snakemake --version` FIRST and branch all execution guidance on 7 vs 8/9.

# Snakemake Workflows

**"Build a reproducible bioinformatics pipeline with Snakemake"** -> Declare each step as a rule that says "a file matching THIS output pattern is produced FROM those inputs", let the engine resolve the DAG backward from requested targets, fan out over samples with wildcards, and pin the software environment so the result reproduces next year.
- Python: Snakefile `rule`/`checkpoint` blocks with `expand()`, `wildcards`, `config`, `resources`, and `conda:`/`container:` (Snakemake)

## The governing principle: Snakemake is pull/goal-oriented - it builds a STATIC DAG backward from requested target files

Snakemake is a pull, make-like engine. An author does NOT describe a forward flow of data. Each `rule` is a pattern-matched recipe ("a file that looks like THIS can be produced FROM that"), and the engine takes the requested target files and works BACKWARD, unifying wildcards by string-matching output filename patterns, until it reaches files already on disk. The whole plan - a static DAG - is computed at parse time, before a single job runs (Köster & Rahmann 2012 *Bioinformatics* 28:2520-2522). Almost every Snakemake bug a biologist hits is a downstream consequence of this one model:

- Rules are wired by OUTPUT-FILE PATTERN, not call order. A missing or typo'd output path silently drops a rule from the DAG - there is no error, the job just never runs. Debugging "why didn't it run" means tracing the backward dependency from the target, not reading top-to-bottom.
- Because the plan is fully known up front, `snakemake -n` (dry run), `--dag`, and `--report` are first-class. This is the payoff of the static model. Nextflow's reactive-dataflow model (processes connected by asynchronous channels, DAG emerges at runtime) has no true dry-run - hold both models in mind and most "why did/didn't it run" questions answer themselves.
- Data-dependent branching is impossible in the base model. If the NUMBER or identity of outputs is unknown until a step runs (split into one file per detected cluster, scatter over however many contigs an assembler emits), the static DAG cannot represent it -> that is exactly what CHECKPOINTS exist for. A biologist who thinks "the pipeline decides at runtime how many chunks" is fighting the paradigm and needs a checkpoint, not a clever `run:` block.
- A workflow manager buys reproducible LOGIC and nothing else automatically. The DAG being deterministic says nothing about tool versions. A rule with no `conda:`/`container:` runs against whatever is on `$PATH`; "reproducible" is unearned until the software environment is pinned (Grüning et al. 2018 *Cell Syst* 6:631-635). Pin containers by DIGEST and conda by LOCKFILE - see Software Deployment below.

## Decision: Snakemake vs Nextflow (pick by team and infrastructure, not benchmarks)

| Dimension | Snakemake | Nextflow | Best when |
|-----------|-----------|----------|-----------|
| Model | pull/make, static DAG at parse time | push/reactive dataflow, dynamic DAG | Snakemake: the plan must be visible before committing an allocation |
| Language | Python DSL (real Python + pandas in the Snakefile) | Groovy DSL | Snakemake: Python-native lab, file-pattern logic |
| Dry run / DAG viz | first-class (`-n`, `--dag`, `--report`) | no true dry-run (`-stub`/`-preview` check wiring only) | Snakemake: HPC where a bad plan is expensive |
| Data-dependent branching | needs checkpoints (escape hatch) | native (channels) | Nextflow: shape depends on runtime data |
| Community pipelines | Workflow Catalog / wrappers (smaller) | nf-core (large, curated) | Nextflow: run a maintained pipeline as-is |
| Sweet spot | single-lab reproducible research, HPC, tight Python integration | cloud/production, multi-institution, nf-core stacks | choose by the ecosystem to integrate with |

Reuse before authoring: for a mainstream analysis (RNA-seq, variant calling, ATAC-seq), a curated community pipeline already encodes years of QC and edge cases. Adopting one means RUNNING it (e.g. nf-core/rnaseq via workflow-management/nf-core-pipelines), not authoring Groovy - so a Python-shop preference for Snakemake only decides the authoring case, not whether to build at all. Author from scratch only for a novel method or an unsupported combination of steps.

## Decision: rerun triggers - why a job reran, or did not

Since Snakemake 7.8 the default is NOT pure mtime. A rerun fires on a SET of triggers: `{mtime, params, input, code, software-env}` (Mölder et al. 2021 *F1000Research* 10:33). This surprises everyone upgrading from old Snakemake.

| Want | Use |
|------|-----|
| classic Make behavior, minimize surprise reruns | `--rerun-triggers mtime` |
| max reproducibility (default) | all five triggers |
| ignore a stable reference's timestamp | `ancient("ref.fa")` on that input |
| mark results current without recompute | `--touch` |
| force specific rules | `--forcerun rule` / `-R` |
| see what WOULD rerun and why | `snakemake -n -R` / `--list-changes code` |

The `code` trigger catches shell/script/run body changes - reformatting whitespace or editing a comment counts as a code change and reruns the job. On very large DAGs or multi-TB inputs the provenance triggers add a hashing/stat storm; `--rerun-triggers mtime` skips it.

## Decision: run vs script vs shell vs notebook vs wrapper

| Situation | Pick | Why |
|-----------|------|-----|
| call a CLI tool (samtools, bwa) | `shell:` | subprocess, conda/container-isolated |
| reusable Python/R analysis needing isolation | `script:` | separate process, `snakemake` object injected |
| standard tool, do not want to write shell | `wrapper:` (PINNED tag) | maintained, ships its own env |
| exploratory, want a re-runnable notebook | `notebook:` | params injected, `--edit-notebook` |
| trivial in-Snakefile glue only | `run:` | NEVER heavy work |

`run:` executes IN the main Snakemake process - it shares the interpreter and GIL, cannot be conda/container-isolated (the `conda:` directive is disallowed with `run:`), blocks the scheduler, and an OOM in it takes down the whole workflow. Move anything beyond trivial glue to `script:`.

## Decision: execution backend (Snakemake 8/9)

| Target | Command |
|--------|---------|
| laptop/workstation | `snakemake --cores N --sdm conda` |
| SLURM (native) | `pip install snakemake-executor-plugin-slurm` then `--executor slurm --jobs N --default-resources` |
| SLURM (legacy sbatch string) | `pip install snakemake-executor-plugin-cluster-generic` then `--executor cluster-generic --cluster-generic-submit-cmd "sbatch ..."` |
| S3/GCS I/O | `pip install snakemake-storage-plugin-s3` then `--default-storage-provider s3 --default-storage-prefix s3://.../` |
| thousands of tiny jobs | add `group:` / `--group-components` to collapse scheduler overhead |

Porting a v7 `--cluster "sbatch --account=X --partition=Y --mem=Z --time=T"` command: for the native `slurm` executor, map those sbatch flags to resource keys (`slurm_account`, `slurm_partition`, `mem_mb`, `runtime` in minutes) set per-rule in `resources:` or globally via `--default-resources`; for a drop-in port keep the old string under `cluster-generic` (its own plugin, above). `--cores` = local cores; `--jobs`/`-j` = number of concurrent cluster/cloud jobs (in v8 these are separate). Profiles are versioned: the file is `config/config.v8+.yaml`, every long option becomes a YAML key.

## Rules, wildcards, and expand

`expand()` returns a LIST of strings by combinatorial substitution - it does NOT touch the filesystem. Use it to enumerate targets in `rule all`.

```python
configfile: 'config/config.yaml'
SAMPLES = config['samples']

rule all:                                              # the requested targets; the DAG is built backward from here
    input:
        expand('results/{sample}.bam', sample=SAMPLES)

rule align:
    input:
        r1 = 'data/{sample}_R1.fq.gz',
        r2 = 'data/{sample}_R2.fq.gz',
        index = 'ref/genome.fa'
    output:
        bam = 'aligned/{sample}.bam'                   # this OUTPUT PATTERN, matched against the target, wires the rule in
    threads: 8
    log:
        'logs/align/{sample}.log'
    shell:
        'bwa mem -t {threads} {input.index} {input.r1} {input.r2} | '
        'samtools sort -@ {threads} -o {output.bam} 2> {log}'
```

## Wildcard constraints - stop silent mis-routing

Wildcards are greedy regex string-unification (`{sample}` compiles to `.+`), not typed parameters. An unconstrained wildcard swallows path separators and adjacent tokens: `data/{sample}.txt` matches `data/a/b.txt` as `sample=a/b`, and `{a}.{b}.txt` on `101.B.normal.txt` has no unique parse. The failure is silent mis-routing, not an error. Constrain whenever a value can contain `/`, `.`, or `_`, or a filename has multiple variable tokens.

```python
wildcard_constraints:
    sample = '[^/]+',                                  # no path separators
    chrom = r'\d+|X|Y|MT'                              # only real chromosome tokens

# two rules whose output patterns can both produce a requested file raise AmbiguousRuleException;
# prefer non-overlapping constraints to disambiguate, and fall back to `ruleorder: a > b` only if needed.
```

## Checkpoints - the ONLY data-dependent-DAG mechanism

**Goal:** produce downstream jobs for a set of files whose number and identity are unknown until a step runs (split a FASTA into one file per detected cluster; scatter over an assembler's contigs).

**Approach:** declare the producing step a `checkpoint` with a `directory()` output; in an input function on the AGGREGATING rule, call `checkpoints.<name>.get(**wildcards)` FIRST - its exception is what forces the engine to run the checkpoint and RE-EVALUATE the DAG - then `glob_wildcards` the checkpoint's declared output dir and `expand()` the real targets.

```python
checkpoint split_fasta:
    input:
        'data/all.fasta'
    output:
        directory('split/{sample}')                    # directory() because the file set is unknowable at parse time
    shell:
        'split_by_cluster.py {input} split/{wildcards.sample}'

def gather_clusters(wildcards):
    # .get() RAISES until the checkpoint has run; that exception drives DAG re-evaluation.
    # Omitting it globs at parse time (empty), so the aggregation silently gets zero inputs - the classic bug.
    ckpt_dir = checkpoints.split_fasta.get(**wildcards).output[0]
    ids = glob_wildcards(f'{ckpt_dir}/{{id}}.fasta').id
    return expand('processed/{sample}/{id}.done', sample=wildcards.sample, id=ids)

rule aggregate:                                         # the input function MUST be attached to the rule that consumes the set
    input:
        gather_clusters
    output:
        'results/{sample}_summary.txt'
    shell:
        'cat {input} > {output}'
```

Point `glob_wildcards` at the checkpoint's declared `directory()` output (a fresh dir) so stale files do not leak into the glob. Prefer one scatter->gather to chains of nested checkpoints.

## Resources, escalating retries, and grouping

Resource callables differ by directive: `resources` is `callable(wildcards [, input] [, threads] [, attempt])`; `threads` is `callable(wildcards [, input])` only; getting the signature wrong is a top error source. `runtime` is in MINUTES. `attempt` starts at 1 and increments per retry - the canonical fix for OOM-killed jobs.

```python
rule call_variants:
    input:
        bam = 'aligned/{sample}.bam'
    output:
        'results/{sample}.vcf'
    threads: 4
    retries: 3                                          # or global --retries 3
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * attempt,   # 8 GB, doubling to 16/24 on OOM restart
        runtime = 240                                   # MINUTES, not seconds; SLURM wall-time
    log:
        'logs/call/{sample}.log'
    shell:
        'variant_caller --threads {threads} {input.bam} > {output} 2> {log}'
```

For thousands of tiny jobs on HPC, per-job scheduler latency dominates: assign rules a `group:` (or `--group-components rule=N`) so they submit as one job. `temp('x.bam')` deletes an intermediate once all consumers are done (huge for disk); `ancient('ref.fa')` excludes an input from mtime-based rerun decisions.

## Software deployment - the layer that makes it reproducible

A clean DAG over unpinned tools is not reproducible. The engine gives layer 1 (logic); the author must pin the software environment. Declare `conda:` (a pinnable file) or `container:` per rule, and activate deployment at run time.

```python
rule fastqc:
    input:
        'data/{sample}.fq.gz'
    output:
        'qc/{sample}_fastqc.html'
    conda:
        'envs/qc.yaml'                                  # a FILE (pinnable), not a bare named env
    container:
        # PIN BY DIGEST, never a mutable tag - :latest or a re-pushed :0.7.17 silently changes the tool and busts the cache
        'docker://quay.io/biocontainers/fastqc@sha256:<digest>'
    shell:
        'fastqc {input} -o qc/'
```

```bash
snakemake --sdm conda --cores 8                         # build per-rule conda envs (was --use-conda in v7)
snakemake --sdm apptainer --cores 8                     # run each rule in its container (was --use-singularity)
snakemake --sdm conda apptainer --cores 8               # containerized conda: build the env INSIDE the pinned image
```

A bare `environment.yml` with `samtools` (no version) resolves differently over time; pin exact builds with a lockfile (conda-lock) for bit-reproducibility. `--containerize` auto-generates a Dockerfile baking all conda envs into one image. Between-workflow caching (`cache: True` + `SNAKEMAKE_OUTPUT_CACHE`) reuses results across workflows but ONLY for deterministic rules - a nondeterministic tool poisons the shared cache with wrong results silently.

## Modularization and reuse

`include: 'rules/x.smk'` is textual inclusion sharing one namespace. For genuine composition of published workflows, use the module system (`module other: snakefile: '...'; use rule * from other as other_*`), which can import, prefix, and override rules. `wrapper: 'v5.0.2/bio/bwa/mem'` pulls a maintained, conda-shipping wrapper - PIN the leading version tag; an unpinned wrapper drifts silently.

```python
include: 'rules/qc.smk'
include: 'rules/align.smk'

rule all:
    input:
        rules.qc_all.input,
        rules.call_all.input
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| A rule silently never runs | its output pattern does not match any requested target (typo/path); the DAG dropped it | trace backward from `rule all`; run `snakemake -n` and inspect the DAG; fix the output path |
| Wildcard captures too much / wrong sample | unconstrained greedy `.+` swallowed a delimiter or path separator | add `wildcard_constraints` (e.g. `sample='[^/]+'`) |
| Aggregation after a split has zero inputs | forgot `checkpoints.X.get()`, or globbed at parse time, or input function on the wrong rule | call `.get(**wildcards)` first, glob the checkpoint's `directory()` output, attach the function to the consumer |
| Everything reruns after a cosmetic edit | the `code` trigger - editing the shell/script body (even whitespace) counts | expected under provenance triggers; use `--rerun-triggers mtime` to opt out |
| Expected a rerun, got none | old mtime mental model, or output is newer than input | check triggers; `--forcerun rule` / `-R` |
| `--cluster`/`S3RemoteProvider` errors on v8 | removed from core in Snakemake 8 | install the executor/storage plugin; use `--executor slurm` and `storage.s3(...)` |
| OOM-killed job (exit 137) | static `mem_mb` too low for the largest sample | `retries` + `mem_mb=lambda wildcards, attempt: base*attempt` |
| Per-job scheduler meltdown on HPC | thousands of tiny jobs, submission overhead dominates | `group:` / `--group-components` to batch into one submission |
| `MissingOutputException` on NFS/Lustre | networked filesystem lags after a job finishes | raise `--latency-wait` |
| Workflow refuses to continue after a killed job ("Incomplete files") | a job died mid-write (SLURM kill, node crash), so its outputs are flagged incomplete | re-run with `--rerun-incomplete` (`--ri`); this is Snakemake's crash-resume, distinct from the rerun triggers |
| Heavy `run:` block hangs the workflow | runs in the main process, shares the GIL, no isolation | move to `script:` |
| "Reproducible" but results differ on a colleague's cluster | no `--sdm`, or a mutable `:latest` container tag | declare `conda:`/`container:`, pin by digest + conda lockfile |

## Related Skills

- workflow-management/nextflow-pipelines - The reactive-dataflow alternative; author here when the DAG must emerge from runtime data
- workflow-management/nf-core-pipelines - Run a curated community Nextflow pipeline instead of authoring from scratch
- workflows/rnaseq-to-de - End-to-end RNA-seq-to-differential-expression pipeline this engine can orchestrate
- read-qc/quality-reports - The QC step a pipeline wraps as an early rule
- read-alignment/bwa-alignment - The alignment tool invoked inside a Snakemake `shell:` rule

## References

- Köster J, Rahmann S. 2012. Snakemake - a scalable bioinformatics workflow engine. *Bioinformatics* 28(19):2520-2522.
- Mölder F, Jablonski KP, Letcher B, et al. 2021. Sustainable data analysis with Snakemake. *F1000Research* 10:33.
- Grüning B, Chilton J, Köster J, et al. 2018. Practical computational reproducibility in the life sciences. *Cell Systems* 6(6):631-635.
- Di Tommaso P, Chatzou M, Floden EW, et al. 2017. Nextflow enables reproducible computational workflows. *Nature Biotechnology* 35(4):316-319.
