---
name: bio-workflow-management-nf-core-pipelines
description: Runs and configures curated nf-core community Nextflow pipelines (rnaseq, sarek, atacseq, methylseq, ampliseq, taxprofiler, fetchngs) reproducibly, pinning the pipeline revision with -r and selecting a container engine and institutional config via -profile. Use when deciding to adopt a community pipeline versus author one from scratch; picking a pipeline and pinning its -r revision; selecting -profile test/docker/singularity/conda plus an institutional config from nf-core/configs; building and validating a samplesheet CSV against the pipeline schema (nf-schema); choosing --genome/iGenomes versus custom references; configuring resources and max_memory for SLURM/AWS Batch; using -resume and -stub; and reading MultiQC outputs.
tool_type: cli
primary_tool: nf-core
---

## Version Compatibility

Reference examples tested with: Nextflow 24.04+, nf-core/tools 3.0+, Docker 24+ or Singularity 3.8+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: `-r <tag>` pins the pipeline to an immutable release; without it `nextflow run nf-core/<pipe>` pulls whatever the mutable default branch is today, so results are not reproducible. Pin the container engine too (a pipeline release ships digest-pinned images; `-profile docker` uses them, `:latest` does not). Schema validation moved from the deprecated nf-validation plugin to nf-schema; a current pipeline uses nf-schema, so validate samplesheets against the pipeline's shipped schema, not a hand-written one.

# nf-core Pipelines

**"Run a curated community pipeline on my samples"** -> Select a versioned nf-core pipeline, pin its release, choose a container profile, validate a samplesheet against the pipeline schema, point it at references, and run it, reading the aggregated MultiQC report at the end.
- CLI: `nextflow run nf-core/<pipeline> -r <version> -profile <container>,<institution> --input samplesheet.csv --outdir results -resume`
- CLI: `nf-core pipelines list` / `nf-core pipelines download` (browse and cache pipelines)

## The governing principle: ADOPT a community pipeline before authoring a new one

For any mainstream analysis - RNA-seq, germline/somatic variant calling, ATAC-seq, ChIP-seq, methylation, amplicon/metagenomics, single-cell - a curated nf-core pipeline already encodes years of QC, edge-case handling, CI/nf-test regression tests, institutional configs for hundreds of HPCs, a standardized samplesheet+schema, and MultiQC reporting (Ewels 2020 *Nat Biotechnol* 38:276-278). Reinventing that in hand-written Nextflow is months of work and worse QC: the community pipeline has already found the bugs a bespoke version will rediscover, and it is maintained across every future tool update and reference build. The decision every biologist should default toward is ADOPT, not BUILD.

The corollary trap is treating "adopt" as "run once and trust the number". A community pipeline is only reproducible if the RUN is pinned: `-r` pins the pipeline version, the release's digest-pinned containers pin the software, `--genome`/reference URIs pin the reference data (Wratten et al. 2021 *Nat Methods* 18:1161-1168; Grüning et al. 2018 *Cell Syst* 6:631-635). An unpinned `nextflow run nf-core/rnaseq` on the default branch with a `:latest` engine is exactly as irreproducible as a hand-rolled script - the curation buys nothing if the invocation is loose. Author from scratch only for genuinely novel logic with no community pipeline, an unsupported combination of steps, or an institutional constraint no config can express (see workflow-management/nextflow-pipelines).

## Decision: adopt an nf-core pipeline vs author a new one

| Situation | Verdict | Why |
|-----------|---------|-----|
| Mainstream analysis with an existing nf-core pipeline (rnaseq, sarek, atacseq, ...) | ADOPT | Curated, CI-tested, institutional configs, MultiQC; DIY is worse QC |
| A supported pipeline plus a few extra params/references | ADOPT + configure | `-c custom.config`, `-params-file`, `--genome`; no authoring needed |
| Genuinely novel method, no community pipeline exists | BUILD | Author in Nextflow; still install tested nf-core modules, do not hand-write wrappers |
| An unsupported ORDER/combination of otherwise-standard steps | BUILD (or fork) | Scaffold with `nf-core pipelines create`; reuse `nf-core modules install` |
| One-off, few linear steps, single sample, single machine | Neither | A plain script is honest; a workflow manager is overhead below this threshold |

## Decision: which container profile by platform

| Platform | Profile | Why |
|----------|---------|-----|
| Laptop / workstation with Docker | `-profile docker` | Simplest; needs root/daemon; digest-pinned images from the release |
| Shared HPC (no root, has Singularity/Apptainer) | `-profile singularity` | Rootless; the default on most academic clusters |
| Cluster allowing Podman | `-profile podman` | Rootless Docker-compatible alternative |
| No container engine available at all | `-profile conda` | Last resort; slower, less reproducible than a pinned image |
| Named institution in nf-core/configs (uppmax, crick, ...) | `-profile singularity,<institution>` | Institutional config sets executor, queues, `max_memory`; comma, no space |
| Smoke test before real data | `-profile test,docker` | Ships a tiny public dataset; proves the install end-to-end in minutes |

Profiles are comma-separated with NO spaces and applied left-to-right (later overrides earlier), so `-profile test,docker` runs the test dataset under Docker, and `-profile singularity,uppmax` layers the institutional config over Singularity.

## Decision: --genome/iGenomes vs custom references

| Reference source | Use when | Caveat |
|------------------|----------|--------|
| `--genome GRCh38` (iGenomes) | A standard build suffices and convenience matters | iGenomes builds are frozen/aging; the annotation may lag current releases |
| Explicit `--fasta` + `--gtf` (+ `--gff`) | A specific build/patch or a non-model organism is needed | Pin the exact reference version; record its URI for provenance |
| Pipeline builds its own index vs `--<tool>_index` | Reusing an index across runs saves hours | A stale index built from a different FASTA silently corrupts results |

The reference layer is a reproducibility layer in its own right: `--genome GRCh38` without a recorded iGenomes snapshot pins less than an explicit `--fasta`/`--gtf` URI pair. Prefer explicit references and record their source when the result must be reproduced.

## The run pattern (pin everything)

```bash
# Smoke test first: tiny public dataset proves the install + engine end to end.
nextflow run nf-core/rnaseq -r 3.14.0 -profile test,docker --outdir results_test

# Real run: -r pins the release, -profile picks the engine, --input is the samplesheet.
nextflow run nf-core/rnaseq -r 3.14.0 \
    -profile singularity \
    --input samplesheet.csv \
    --genome GRCh38 \
    --outdir results \
    -resume
```

- `-r 3.14.0` is the pipeline REVISION (a git tag). It is MANDATORY: without it the run tracks the mutable default branch and is not reproducible.
- `-profile singularity` selects the container engine (comma-add an institutional config: `-profile singularity,uppmax`).
- Single-dash options (`-r`, `-profile`, `-resume`, `-c`, `-params-file`) are NEXTFLOW options; double-dash options (`--input`, `--genome`, `--outdir`, `--max_memory`) are PIPELINE parameters. Mixing up the dash count is the most common invocation error.
- Config precedence, low to high: the pipeline's built-in `nextflow.config` -> `conf/base.config` -> selected profiles -> `-c custom.config` -> `-params-file params.yaml` -> a `--param` on the command line. A later source overrides an earlier one.

Parameters can be supplied in a YAML/JSON file instead of long command lines, which is the reproducible-provenance form:

```yaml
# params.yaml  (nextflow run nf-core/rnaseq -r 3.14.0 -profile singularity -params-file params.yaml)
input: samplesheet.csv
genome: GRCh38
outdir: results
aligner: star_salmon
```

## The samplesheet and schema validation

Every nf-core pipeline reads a CSV samplesheet whose exact columns are defined by the pipeline's shipped schema (`assets/schema_input.json`) and validated by the nf-schema plugin at launch, before any compute is spent. Wrong or misordered columns fail fast with a schema error rather than hours in.

```csv
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,/data/ctrl1_R1.fastq.gz,/data/ctrl1_R2.fastq.gz,auto
CONTROL_REP2,/data/ctrl2_R1.fastq.gz,/data/ctrl2_R2.fastq.gz,auto
TREAT_REP1,/data/treat1_R1.fastq.gz,/data/treat1_R2.fastq.gz,auto
```

- Column names are pipeline-specific: nf-core/rnaseq uses `sample,fastq_1,fastq_2,strandedness`; nf-core/sarek uses `patient,sample,lane,fastq_1,fastq_2`. Read the pipeline's `docs/usage.md` for the exact schema, never guess.
- A single-end sample leaves `fastq_2` empty; multiple rows sharing one `sample` value are merged (technical replicates / multiple lanes), which is how the pipeline knows to concatenate them.
- Absolute paths or URLs are safest; relative paths resolve against the launch directory.

Threading columns into the pipeline is handled by the META MAP convention (worth understanding when reading logs or outputs): each sample flows internally as a tuple `[ meta, files ]` where `meta` is a map like `[ id:'CONTROL_REP1', single_end:false ]`. The samplesheet columns become `meta` keys, so sample identity and pairing travel WITH the files through every step - which is why outputs and the MultiQC report are labelled by the `sample` value from the sheet.

## Institutional configs and resource limits

nf-core/configs supplies ready-made profiles for hundreds of clusters (executor, queues, module system, resource ceilings). Use a named one when it exists; otherwise write a small custom config.

```groovy
// custom.config  (nextflow run ... -c custom.config)
process {
    executor = 'slurm'
    queue    = 'normal'
    // Clamp per-process resource escalation to the real node ceiling, so an auto-retry that
    // doubles memory never requests more than a node has. resourceLimits is the nf-core/tools 3.0+
    // form (Nextflow 24.04+); it replaced the deprecated params.max_cpus/max_memory/max_time + check_max().
    resourceLimits = [ cpus: 32, memory: 128.GB, time: 48.h ]
}
```

- nf-core pipelines escalate resources on retry (a task that OOM-kills retries with more memory); `process.resourceLimits` caps that escalation so a request stays schedulable. A pipeline built on the pre-3.0 template instead reads `params.max_cpus`/`max_memory`/`max_time` (the `check_max()` pattern, deprecated and removed from the template in tools 3.0) - match whichever the pinned `-r` release ships. Set the ceiling to the real node/queue limits either way.
- Do not edit the pipeline's own `conf/base.config`; layer overrides through `-c custom.config` so the pipeline stays a clean, updatable checkout.

## Resume, stub, and previewing the plan

```bash
# -resume reuses cached tasks whose inputs+script+container hash is unchanged.
nextflow run nf-core/rnaseq -r 3.14.0 -profile singularity --input samplesheet.csv --outdir results -resume

# -stub runs each process's stub block (touch fake outputs) to validate wiring in seconds.
nextflow run nf-core/rnaseq -r 3.14.0 -profile test,docker --outdir results -stub
```

- `-resume` keys on a hash of each task's inputs, resolved script, and container reference. On a network filesystem (Lustre/NFS) unreliable mtimes cause spurious cache misses; the standard fix is `cache 'lenient'` in a custom config. Deleting `work/` destroys the resume cache - a re-run then recomputes everything.
- `-stub` validates that the samplesheet, profile, and channel wiring are correct without running any tool, which is the fast pre-flight before committing an HPC allocation.

## Building a new pipeline (only when adoption does not fit)

```bash
# Scaffold a standardized pipeline (template, CI, lint, nf-test) - current tools syntax.
nf-core pipelines create

# Install a pre-written, tested module instead of hand-writing a tool wrapper.
nf-core modules install fastqc
nf-core subworkflows install bam_sort_stats_samtools

# Lint against the template and run the module's nf-test snapshot tests.
nf-core pipelines lint
nf-core modules test fastqc
```

Even when building, reuse the community's tested modules rather than hand-writing bwa/samtools/fastqc wrappers. Authoring mechanics (channels, DSL2, resume internals) live in workflow-management/nextflow-pipelines.

## Interpreting MultiQC output

Every nf-core run aggregates per-tool QC into a single `multiqc_report.html` under the output directory (plus parsed `multiqc_data/` tables). Read it before trusting any downstream result:

- The General Statistics table is per-sample; scan for an outlier column (low aligned %, high duplication, skewed GC, adapter content) that flags a failed library BEFORE it contaminates differential analysis.
- Section order mirrors the pipeline steps (e.g. FastQC -> trimming -> alignment -> quantification for rnaseq); a section missing for one sample means that sample failed a step - cross-check the Nextflow log.
- MultiQC reports what the tools measured; it does not decide pass/fail. Set thresholds from the assay, and treat the report as the triage surface, not the verdict (read-qc/quality-reports).

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Results differ between runs / cannot reproduce a published run | no `-r`, so the mutable default branch was used | always pin `-r <version>`; record it alongside results |
| `Unknown configuration profile` or only one profile applied | `-profile test docker` with a space | use a comma, no space: `-profile test,docker` |
| Launch fails immediately with a schema/validation error | samplesheet columns wrong, misordered, or misnamed for this pipeline | match the pipeline's `assets/schema_input.json` / `docs/usage.md` exactly |
| `--input` or `--genome` "is not a valid parameter" | used a single dash (`-input`) - that is a Nextflow option namespace | pipeline params take double dash; Nextflow options (`-r`, `-profile`, `-resume`) take single |
| `-resume` re-runs everything on the cluster | mtime-based cache misses on a network filesystem | add `cache 'lenient'` via `-c custom.config`; never delete `work/` |
| Container/tool "command not found" at runtime | no container engine profile selected (bare `nextflow run`) | add `-profile docker`/`singularity`/`conda` |
| Task unschedulable, requests more memory than any node | retry escalation exceeded the node ceiling | set `process.resourceLimits = [cpus:, memory:, time:]` (pre-3.0 pipelines: `params.max_memory`/`max_cpus`/`max_time`) |
| Wrong/aging annotation with `--genome` | iGenomes builds are frozen and can lag current releases | supply explicit `--fasta`/`--gtf` for a specific build and record the URI |

## Related Skills

- workflow-management/nextflow-pipelines - Author a Nextflow pipeline from scratch when no community pipeline fits
- workflow-management/snakemake-workflows - Rule-based alternative engine for pipeline authoring
- workflows/rnaseq-to-de - Take an nf-core/rnaseq count matrix into differential expression
- read-qc/quality-reports - Interpret the FastQC/MultiQC QC surface a pipeline emits

## References

- Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. 2020. The nf-core framework for community-curated bioinformatics pipelines. *Nat Biotechnol* 38(3):276-278.
- Di Tommaso P, Chatzou M, Floden EW, Prieto Barja P, Palumbo E, Notredame C. 2017. Nextflow enables reproducible computational workflows. *Nat Biotechnol* 35(4):316-319.
- Wratten L, Wilm A, Göke J. 2021. Reproducible, scalable, and shareable analysis pipelines with bioinformatics workflow managers. *Nat Methods* 18:1161-1168.
- Grüning B, Chilton J, Köster J, et al. 2018. Practical computational reproducibility in the life sciences. *Cell Syst* 6(6):631-635.
