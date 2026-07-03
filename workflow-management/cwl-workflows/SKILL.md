---
name: bio-workflow-management-cwl-workflows
description: Authors portable, strongly-typed bioinformatics pipelines in the Common Workflow Language (CWL v1.2) as CommandLineTool/Workflow/ExpressionTool documents, validated with cwltool and run at scale on Toil/Arvados/Calrissian. Use when deciding CWL (portability/provenance/regulated) vs Nextflow/WDL/Snakemake; declaring secondaryFiles for indexed companions (.bai/.fai/.dict/.tbi and the caret rule); putting resources/containers under requirements (must-hold) vs hints (advisory) to avoid silent OOM; choosing scatterMethod (dotproduct vs flat_/nested_crossproduct); preferring $(...) parameter refs over ${...} JavaScript for portability; pinning DockerRequirement images; or emitting a CWLProv provenance object for audited/clinical settings.
tool_type: cli
primary_tool: cwltool
---

## Version Compatibility

Reference examples tested with: cwltool 3.1+, CWL spec v1.2, Docker 24+ (or Singularity/Apptainer 3.8+)

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: target `cwlVersion: v1.2` for genomics (the record-form secondaryFiles with `required:` needs it). cwltool is the REFERENCE runner - correct, local, single-node, deliberately slow - NOT the spec; run Toil/Arvados/Calrissian at scale. In v1.2, ExpressionTool outputs are NOT type-checked (a known reference-impl gap, fix planned for v1.3), so do not lean on ExpressionTool for type safety.

# CWL Workflows

**"Write a portable pipeline that runs the same everywhere and proves what it ran"** -> Describe each tool and their wiring as strongly-typed, machine-checkable CWL documents that any conforming runner honors identically, validate the contract before any compute, then execute locally (cwltool) or at scale (Toil/Arvados/Calrissian).
- CLI: `cwltool --validate wf.cwl`, `cwltool wf.cwl job.yml`, `cwltool --provenance ro/ wf.cwl job.yml`
- YAML: `class: CommandLineTool` (wrap a tool) and `class: Workflow` (wire tools) at `cwlVersion: v1.2`

## The governing principle: CWL is a SPECIFICATION, not an engine

CWL deliberately splits the workflow DESCRIPTION (a portable, declarative, strongly-typed YAML/JSON document) from its EXECUTION (performed by any conforming runner). Nextflow and Snakemake are engines that happen to have a DSL; WDL is a language with a dominant engine (Cromwell). CWL alone is a community-governed open standard with multiple independent implementations and a formal conformance test suite (Crusoe 2022 *Commun ACM* 65(6):54-63). The whole value proposition - portability, auditability, vendor-neutrality, provenance, regulatory fit - is a consequence of "spec not engine." The correct frame: the author writes a portable, machine-checkable CONTRACT for a computation that any conforming platform must honor identically. That contract is why CWL is the most VERBOSE and most EXPLICIT of the four systems - the verbosity buys static analyzability and portability.

The most common conceptual error is conflating CWL with cwltool. cwltool is the reference implementation: correct, single-node, and slow by design (it prioritizes spec-conformance over speed). "CWL is slow" or "CWL can't scale" almost always means "I ran cwltool" - at scale the same unchanged document runs on Toil (HPC/cloud batch), Arvados (clinical/enterprise data management), or Calrissian (Kubernetes). Never equate a runner's limits with the spec's.

Adopting CWL buys reproducible workflow LOGIC and nothing else automatically. A clean typed DAG over unpinned tools is NOT reproducible - pin the software environment (containers by digest, not a moving `:latest`), the reference data and seeds, and control arch/thread/locale leaks separately. The type system closes the wiring-error layer; the author still owns the rest.

## Decision: choose CWL vs the other engines

| Dimension | CWL | Nextflow | WDL | Snakemake |
|-----------|-----|----------|-----|-----------|
| Nature | open SPEC, many engines | engine + Groovy DSL | language + Cromwell/miniwdl | engine + Python DSL |
| Typing | strong static (File/Dir/record/enum/optional) | dynamic | moderate | weak (paths/strings) |
| Index companions | `secondaryFiles` = File type property | manual channel/tuple wiring | manual per-input | manual |
| Provenance | CWLProv RO out of the box | report/trace/DAG | via platform | report/plugins |
| Verbosity | HIGHEST (deliberate) | terse | moderate | moderate |
| New-author mindshare | DECLINING | ASCENDANT (nf-core) | strong on Terra | strong in academia |
| Sweet spot | multi-platform, provenance-critical, regulated/clinical, standards-driven | fast authoring, curated catalog | Terra/GATK ecosystem | single-lab Python/HPC |

Choose CWL when vendor-neutral portability across multiple platforms, a standardized provenance artifact for a regulated/audited setting, static type-checking before compute, or publishing to a multi-runner registry (Dockstore) drives the decision. Be honest about mindshare: for NEW pipeline authoring CWL has been losing ground to Nextflow/nf-core and WDL/Terra for years (they win on authoring speed and community momentum). CWL retains and deepens its hold where the SPEC is the point - not as a neutral default. Prefer Nextflow/WDL when authoring speed, an existing curated catalog (nf-core), or a specific hosted platform (Terra) dominates.

## Decision: which runner

| Engine | For | Runtime model | Notes |
|--------|-----|---------------|-------|
| cwltool | authoring, `--validate`, `--pack`, `--provenance`, CI, local dev | local, single-node | the conformance yardstick; slow by design - do NOT read its limits as CWL's |
| Toil | HPC and cloud batch scale | Python; Slurm/Kubernetes/AWS/Grid Engine | `toil-cwl-runner`; the workhorse for large CWL |
| Arvados | enterprise/clinical data management + execution | cluster + content-addressed storage | `arvados-cwl-runner`; strong data provenance; regulated settings |
| Calrissian | CWL on Kubernetes | one k8s pod per step | needs `ReadWriteMany` volumes; cloud-native parallelism |
| Cromwell | primarily WDL, PARTIAL CWL | JVM | runs a subset only; do not rely on it for full CWL conformance |

## Decision: requirements vs hints (get this wrong and jobs silently OOM)

`requirements` MUST be satisfied - if the runner cannot honor one, execution FAILS loudly (correctly). `hints` are advisory: the runner MAY honor or ignore them without error. The canonical failure is putting `ResourceRequirement: {ramMin: 32000}` under `hints` on a memory-hungry step - a runner is free to ignore a hint and schedule it on a small node, giving intermittent OOM kills. Anything whose absence would corrupt results or crash (the container, minimum RAM/cores, a required input layout, an env var a tool depends on) goes under `requirements`. Both inherit Workflow -> step -> tool with the INNERMOST declaration winning, so set a default `DockerRequirement`/`ResourceRequirement` at workflow scope and override per-step where a tool needs a different image or more RAM. A surprising container at a step is usually a forgotten override.

## Decision: scatterMethod

`scatter` runs a step once per array element; when scattering over MULTIPLE inputs, `scatterMethod` decides how they combine. This is the most misunderstood CWL construct.

| scatterMethod | Combines by | Jobs | Output shape | Use when |
|---------------|-------------|------|--------------|----------|
| `dotproduct` | position-aligned zip | N (arrays MUST be equal length) | flat array of N | paired arrays that correspond 1:1 (R1[i] with R2[i]) |
| `flat_crossproduct` | every combination | N x M | FLAT array of N x M | all pairs, want a flat result list |
| `nested_crossproduct` | every combination | N x M | NESTED array (N of M) | all pairs, preserve the 2-D grid |

`dotproduct` requires equal-length arrays (unequal is an error, not truncation). The two cross-products run the SAME N x M jobs and differ only in output nesting - choosing `flat_` vs `nested_` wrong gives the right computations with a mis-shaped output that then mis-wires a downstream `File[]` step. `ScatterFeatureRequirement` must be declared regardless of scatterMethod.

## CommandLineTool: wrap one tool

A CommandLineTool binds typed inputs to the command line (`inputBinding`) and captures outputs (`outputBinding.glob` or `stdout`). Tool outputs use `outputBinding`; workflow outputs use `outputSource` - mixing them is a validation error.

```yaml
cwlVersion: v1.2
class: CommandLineTool
baseCommand: [bwa, mem]
requirements:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bwa:0.7.17--he4a0461_11   # digest-pin in production; :latest breaks reproducibility
  ResourceRequirement:            # under requirements: must hold, or the runner fails loudly (a hint could be ignored -> OOM)
    coresMin: 8
    ramMin: 16000                 # MB; bwa-mem index residency + reads, empirical floor for a human genome
inputs:
  reference:
    type: File
    secondaryFiles:               # the .amb/.ann/.bwt/.pac/.sa BWA index files are STAGED next to the fasta
      [.amb, .ann, .bwt, .pac, .sa]
    inputBinding: {position: 2}
  reads_1: {type: File, inputBinding: {position: 3}}
  reads_2: {type: File?, inputBinding: {position: 4}}   # File? is optional: [null, File]
  threads: {type: int, default: 8, inputBinding: {prefix: -t, position: 1}}
stdout: aligned.sam
outputs:
  sam: {type: stdout}
```

## Workflow: wire tools by explicit typed connections

CWL is PULL/goal-oriented and fully declarative: dataflow is wired explicitly through `source`/`outputSource`, never through implicit channels (Nextflow) or filename wildcards (Snakemake). `--validate` type-checks every connection statically, before a byte of data moves.

```yaml
cwlVersion: v1.2
class: Workflow
requirements:
  ScatterFeatureRequirement: {}
inputs:
  fastq_1: File
  fastq_2: File
  salmon_index: Directory
outputs:
  quant_results:
    type: Directory
    outputSource: salmon/quant_dir     # workflow output wires with outputSource (NOT outputBinding)
steps:
  fastp:
    run: fastp.cwl
    in: {reads_1: fastq_1, reads_2: fastq_2}
    out: [trimmed_1, trimmed_2, json_report]
  salmon:
    run: salmon_quant.cwl
    in: {index: salmon_index, reads_1: fastp/trimmed_1, reads_2: fastp/trimmed_2}   # source: other_step/output
    out: [quant_dir]
```

## secondaryFiles: index and companion files as a type property

Genomics tools demand companion files that must sit next to the primary with a derived name: `.bam` needs `.bai`, `.fasta` needs `.fai` and `.dict`, `.vcf.gz` needs `.tbi`. CWL makes the companion a PROPERTY of the File type, so every conforming runner is OBLIGATED to co-stage them - the other engines leave this to hand-wiring. This is the single most genomics-relevant CWL feature and the strongest reason to target v1.2.

```yaml
inputs:
  bam:
    type: File
    secondaryFiles: [.bai]         # append: sample.bam -> sample.bam.bai staged alongside
  reference:
    type: File
    secondaryFiles:                # v1.2 record form makes required-ness explicit
      - {pattern: .fai, required: true}
      - {pattern: ^.dict, required: false}   # caret ^ STRIPS one extension: genome.fasta -> genome.dict (NOT genome.fasta.dict)
```

The caret `^` removes one extension from the basename before appending; each leading `^` strips one more - the classic `.dict` gotcha. In v1.0/v1.1 a bare string pattern is required-by-default; only the v1.2 record `{pattern, required}` form can mark an index optional (a `.tbi` that may be absent). secondaryFiles are declared on OUTPUT File parameters too, so a produced BAM carries its `.bai` to the next step; on inputs required defaults true, on outputs the index is collected if present.

## Expressions: $(...) is portable, ${...} is a portability debt

Parameter references `$(...)` are a restricted, safe subset - property/index access into `inputs`, `self`, `runtime` (e.g. `$(inputs.reads.nameroot)`, `$(runtime.outdir)`). They need NO JavaScript engine, are statically analyzable, and are portable - prefer them. Full JavaScript `${ return ...; }` runs only when `InlineJavascriptRequirement` is present; it is powerful but (a) needs a node engine wherever the workflow runs, (b) is opaque to static analysis, (c) is the leading cause of "works on my runner, breaks on theirs." Reach for `${...}` only when a parameter reference genuinely cannot express the need, and treat every `InlineJavascriptRequirement` as a portability debt taken on knowingly. `valueFrom` transforms a step input before it reaches the tool (needs `StepInputExpressionRequirement` for expressions at step scope).

## Portability leaks: "run anywhere" is real but DISCIPLINED

The promise is genuine for disciplined CWL, but it leaks - name the leaks. `InlineJavascriptRequirement` needs a JS engine (node). `DockerRequirement` carries arch/registry assumptions (an amd64-only image fails on arm64 Apple Silicon/Graviton; a private-registry image needs credentials the target may lack; `:latest` is not reproducible). Engine-specific extensions under custom namespaces (`cwltool:`, `arv:`) are portable only among runners that understand them - they live under `hints` so a naive runner can ignore them, but a workflow that DEPENDS on their behavior has forfeited portability. Conformance is graded (a coverage percentage per implementation), not binary - "valid CWL" does not guarantee "runs identically on engine X." The disciplined recipe: target v1.2; containerize every tool with a digest-pinned multi-arch image; minimize `${...}` in favor of `$(...)`; keep `cwltool:`/`arv:` items under `hints` and never depend on them for correctness; in the input (job) object prefer `location:` URIs over a local `path:` (a `path:` binds the job to one machine's filesystem); validate, then exercise the real target engine before trusting portability.

## Provenance: CWL's high ground for regulated/clinical settings

`cwltool --provenance ro/ wf.cwl job.yml` produces a CWLProv Research Object (a W3C PROV + RO-Crate/BagIt bundle) capturing the workflow, the exact input object, all outputs, intermediates, container images, and the enactment trace (Khan 2019 *GigaScience* 8(11):giz095). No other mainstream system ships a standardized retrospective-provenance artifact out of the box - this is the concrete reason CWL wins in regulated/audited genomics: hand an auditor one object that answers "exactly what ran, on what inputs, in what containers, producing what outputs." CWL is also a first-class GA4GH citizen: TRS (Tool Registry Service, implemented by Dockstore) standardizes discovery, WES standardizes cross-platform execution - the strongest reproducibility+portability+provenance story of the four systems.

## Run commands

```bash
cwltool --validate wf.cwl                     # static type-check the contract; no compute
cwltool wf.cwl job.yml                         # run locally (reference runner)
cwltool --singularity wf.cwl job.yml           # swap container runtime (Apptainer is Singularity-compatible)
cwltool --pack wf.cwl > packed.cwl             # bundle a multi-file workflow into one shareable JSON
cwltool --provenance ro/ wf.cwl job.yml        # emit a CWLProv Research Object
toil-cwl-runner --batchSystem slurm wf.cwl job.yml   # same document, at HPC scale
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| "CWL is slow / can't scale" | ran cwltool (reference runner, single-node, slow by design) | run the same document on Toil/Arvados/Calrissian; do not equate runner limits with the spec |
| "index not found" at runtime (e.g. no `.bai`/`.fai`) | secondaryFiles not declared, so the runner staged the primary but not its index | declare `secondaryFiles` on the indexed File input (and output) |
| `genome.fasta.dict` produced instead of `genome.dict` | forgot the caret; `.dict` appends, `^.dict` strips one extension | use `^.dict` (each `^` strips one extension) |
| Intermittent OOM / undersized node on a heavy step | `ResourceRequirement` placed under `hints` (advisory, may be ignored) | move anything that MUST hold under `requirements` |
| Validation error on a workflow output | used `outputBinding.glob` on a workflow output | workflow outputs wire via `outputSource: step/out`; only tool outputs use `outputBinding` |
| Scatter runs but produces N x M or wrong-shaped output | wrong scatterMethod (dotproduct vs flat_/nested_crossproduct) | pick deliberately: dotproduct=zip, flat_/nested_=all pairs (flat vs nested output) |
| Works on cwltool, fails/differs on another engine | `${...}` JS or `cwltool:`/`arv:` extension the target lacks; graded conformance | prefer `$(...)`; keep engine extensions in hints; test the real target |
| "ScatterFeatureRequirement not specified" | used `scatter:` without the feature flag | add `requirements: [ScatterFeatureRequirement]` |
| Non-reproducible result across time | mutable `:latest` container tag | pin `DockerRequirement` by `@sha256:` digest |

## Related Skills

- workflow-management/wdl-workflows - WDL/Cromwell alternative for the Terra/GATK ecosystem
- workflow-management/nextflow-pipelines - reactive-dataflow alternative with the nf-core catalog
- workflow-management/snakemake-workflows - Python/file-pattern alternative for single-lab HPC
- workflows/fastq-to-variants - an end-to-end variant-calling pipeline these engines orchestrate

## References

- Crusoe MR, Abeln S, Iosup A, et al. 2022. Methods Included: Standardizing Computational Reuse and Portability with the Common Workflow Language. *Commun ACM* 65(6):54-63. DOI 10.1145/3486897.
- Amstutz P, Crusoe MR, Tijanic N, et al. 2016. Common Workflow Language, v1.0. figshare. DOI 10.6084/m9.figshare.3115156.v2.
- Khan FZ, Soiland-Reyes S, Sinnott RO, Lonie A, Goble C, Crusoe MR. 2019. Sharing interoperable workflow provenance: a review of best practices and their practical application in CWLProv. *GigaScience* 8(11):giz095. DOI 10.1093/gigascience/giz095.
- Wratten L, Wilm A, Göke J. 2021. Reproducible, scalable, and shareable analysis pipelines with bioinformatics workflow managers. *Nat Methods* 18:1161-1168. DOI 10.1038/s41592-021-01254-9.
