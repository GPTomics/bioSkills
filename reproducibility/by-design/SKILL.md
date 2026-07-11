---
name: bio-reproducibility-by-design
description: Designs bioinformatics workflows for reproducibility and reuse from the start, treating code-on-GitHub-with-a-Docker-image as availability, not reproducibility. Use when deciding which of the three levels of reproducibility (L1 same data+code+environment, L2 robustness, L3 scientific) a claim or reviewer request means; restructuring a monolithic script into modules with explicit input/output contracts, no hard-coded paths, and no data/code mixing so blocks are independently reusable; establishing tool identity across papers and code (bio.tools/EDAM, versioning, CoPaLink entity linking) so two workflows provably use the same tool; reconstructing and visualizing hand-written workflow structure with BioFlow-Insight to keep the diagram synced with the code; diagnosing workflow decay; and deciding whether to adopt a curated community pipeline given that almost nothing bespoke is reused as-is. Routes environment pinning, testing, and provenance to the sibling skills.
tool_type: mixed
primary_tool: BioFlow-Insight
---

## Version Compatibility

Reference examples tested with: BioFlow-Insight 0.2+, Nextflow 24.04+, Snakemake 8+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show bioflow-insight` then `bioflow-insight --help` to check the CLI
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: BioFlow-Insight parses Nextflow DSL2 (not DSL1, removed in 22.12) and is under active
development, so its CLI/output flags move between releases; treat the commands below as the shape,
and confirm against `--help`. bio.tools and EDAM are living registries queried over HTTP.

# Reproducibility by Design

**"Make this analysis reproducible and reusable, not just available"** -> Design the workflow so a stranger can re-run it in a pinned environment on reference data, understand its structure without reading every line, and trust that each output traces to the inputs that produced it. Availability (a Git repo + a Docker image) is a prerequisite, not the goal.

## The governing principle: availability != reproducibility != reuse

Three distinct properties, routinely conflated under editorial and reviewer pressure (Cohen-Boulakia 2017 *Future Gener Comput Syst* 75:284-298):

- **Availability** is that the code and an image exist somewhere. It is necessary and cheap, and where most projects stop. Code dumped this way is hard to read: the analysis steps are hard to guess, several authors mixed styles, data and code are entangled (paths hard-coded mid-script, which itself breaks reproducibility), and the docs record *what* was done but almost never *why* and *how*.
- **Reproducibility** is the ability to regenerate the result. It requires the pinned environment, the reference data and parameters, and enough structure to be understood. "Trust me" is not a scientific stance; understanding requires structure, and raw code does not provide it.
- **Reuse** is others building on the parts. It is impossible without reproducibility: you cannot reuse what you cannot re-run and comprehend. Publishing a new-method-vs-old-method comparison without first checking the old pipeline still runs identically on the old data rests the whole comparison on sand.

The design consequence: reproducibility is not a step added at submission, it is a property built in from the first module. Design for it, or it will not appear later.

## The five pillars (and where each is handled)

| Pillar | What it means | Where |
|--------|---------------|-------|
| 1. Modularity | Step-by-step design, explicit expected input/output per module, blocks reusable independently | this skill |
| 2. Environment & dependencies | Qualify and package each step so it runs on any machine (container by digest, conda by lockfile) | workflow-management (mature tooling exists) |
| 3. Testing | Ship a sample input + expected output; alert automatically when a tool change silently breaks it | reproducibility/workflow-testing (the community's weakest practice) |
| 4. Tool identity & naming | Name, version, and describe each tool's function/inputs so two workflows provably use the same tool | this skill |
| 5. Fine-grained provenance | Link each output only to the inputs that genuinely produced it, without provenance overload | reproducibility/workflow-provenance |

This skill owns pillars 1 and 4 (the design and comprehension layer) and routes the rest. A clean pillar-1 DAG over unpinned tools (no pillar 2) is still not reproducible.

## Pillar 1: modularity in practice

The unit of reuse is the *step* (a "process" in Nextflow, a "rule" in Snakemake), not the pipeline. For each module you must be able to answer three questions: what does it do, how is it used, and what does it depend on. If a consumer cannot answer these without reading the internals, the boundary is wrong.

| Design rule | Reproducibility-by-design rationale | Anti-pattern it kills |
|-------------|-------------------------------------|-----------------------|
| One tool per module, named for its function | The step is understandable and swappable in isolation | 300-line "do everything" script |
| Declare expected input and output explicitly (typed where possible) | The contract is the interface; consumers reason without the internals | guessing I/O from glob patterns |
| No hard-coded paths; inputs via params/channels | The module runs on another machine and another dataset | `read.csv("/home/me/data/x.csv")` mid-script |
| Separate data from code | Data lives in config/samplesheets, not in logic | sample list pasted into the loop body |
| Pin the environment per module (pillar 2) | The step's software identity is fixed | `container: 'tool:latest'` |
| Keep structure and code in one artifact | The diagram cannot drift from what runs (see BioFlow-Insight below) | a slide drawn once, never updated |

An empirical caution to calibrate expectations: across ~10,000 published workflow steps, almost nothing is reused *as-is*; the few reused parts are small glue steps, and even those need adjustment (Cohen-Boulakia 2017). Curated, reviewed workflows are reused substantially more. Two consequences: (1) design a module for reuse deliberately, it does not happen by accident; (2) for any mainstream analysis, adopt a curated community pipeline (nf-core, WARP, WorkflowHub) before authoring, because it already encodes the modularity, tests, and configs (Langer 2025 *Genome Biol* 26:228). See workflow-management for the adopt-vs-author decision.

```python
# modules/dedup.nf equivalent contract, made explicit and reusable.
# GOOD: the interface (in/out) is declared; the tool is pinned; nothing is hard-coded.
process MARK_DUPLICATES {
    tag "${meta.id}"                                    // sample identity threads through
    container 'quay.io/biocontainers/gatk4:4.5.0.0--py36hdfd78af_0'  // pinned, not :latest

    input:
    tuple val(meta), path(bam)                          // declared input contract

    output:
    tuple val(meta), path("*.markdup.bam"), emit: bam   // declared output contract
    path "*.metrics.txt",                    emit: metrics

    script:
    """
    gatk MarkDuplicates -I ${bam} -O ${meta.id}.markdup.bam -M ${meta.id}.metrics.txt
    """
}
```

## Pillar 4: tool identity and naming

Making sure two workflows refer to *the same tool* is remarkably hard, and it is the same class of problem as gene naming solved a generation earlier: without a controlled identity, a name is ambiguous across contexts. A tool called `bwa` in a paper may be invoked as `bwa-mem2` in the code, at an unstated version, wrapped in a renamed container.

| Practice | Mechanism | Payoff |
|----------|-----------|--------|
| Register/resolve tool identity | bio.tools ID; EDAM ontology for operation/format/topic | a stable, machine-resolvable name, not a free-text string |
| Pin and record the version | semantic version + container digest (pillar 2) | the identity cannot drift between paper and re-run |
| Describe function and I/O | EDAM operation + input/output data types | consumers know what the step does without running it |
| Reconcile names across paper and code | CoPaLink: named-entity recognition + entity linking over articles and Nextflow code (Sebe 2026) | the tool a paper *describes* is linked to the tool the code *calls* |

CoPaLink's finding matters for method choice: large language models do not answer the tool-linking task best; simpler, far cheaper NER + entity-linking models outperform them on a corpus of 205 articles and 3,000+ Nextflow workflows. Naming reproducibility is itself a reproducibility problem (Cuevas Villarmin 2024): NER results shift with seed and library version, so pin the extraction pipeline too.

## Workflow comprehension: keep the diagram synced with the code

Hand-written workflow systems (Nextflow, Snakemake) have no built-in graphical view, so people redraw their pipeline by hand in slides. Those drawings look alike but are not: barely-formatted, over-coloured, impossible to keep in sync with evolving code, and untraceable back to the lines that implement them. A stale diagram is worse than none.

BioFlow-Insight (Marchment 2024 *NAR Genom Bioinform* 6(3):lqae092) parses Nextflow code to its abstract syntax tree, reconstructs the workflow structure automatically, and renders an interactive view where clicking a step reveals the underlying code and moves between levels of detail. Because the structure is *derived from the code*, it cannot drift. In a 16-participant study (Master's students and Nextflow developers) on step-counting and step-ordering tasks, users with the tool were faster, more accurate, and more confident, even on tasks they rated hard.

```bash
# Reconstruct and export a workflow's structure straight from the code (shape; confirm flags with --help).
# The point is that the graph is generated, not hand-drawn, so it stays true to what runs.
bioflow-insight main.nf --output-dir structure/         # parse AST -> reconstruct DAG
# Snakemake's native equivalent for the DAG (also code-derived):
snakemake --dag        | dot -Tsvg > dag.svg            # per-job DAG
snakemake --rulegraph  | dot -Tsvg > rulegraph.svg      # rule-level structure (closer to modules)
# Nextflow's native DAG (code-derived):
nextflow run main.nf -preview -with-dag flowchart.mmd   # Mermaid; -preview does not run tasks
```

## Workflow decay: what breaks a re-run over time, and the fix

| Era | Decay cause | Modern state | Design fix |
|-----|-------------|--------------|------------|
| ~20 yr ago | Dependence on external web services that vanished | mostly gone (younger authors never met it) | avoid live-service calls in analysis logic |
| now | Library / container drift when a dependency changes | you update the workflow, which is a better failure mode | pin container by digest, conda by lockfile (pillar 2) |
| now | Re-execution impossible when data cannot be shared | containerization solves the *software* side | ship shareable synthetic reference data (workflow-testing) |
| always | Structure diverges from code; nobody reads the diagram | reconstructable from code | generate the diagram (BioFlow-Insight / --rulegraph), never hand-draw it |

## Decision: what "make it reproducible" should trigger

| The ask | The pillar(s) | Route to |
|---------|---------------|----------|
| "It's on GitHub with Docker, isn't that reproducible?" | 2 only, not 1/3/4/5 | explain the gap; do the rest here |
| "Others can't reuse my code" | 1 (modularity), 4 (naming) | this skill |
| "I can't see what my pipeline does" | comprehension | BioFlow-Insight / rulegraph, this skill |
| "Which tool/version did the paper actually use?" | 4 | bio.tools/EDAM + CoPaLink, this skill |
| "Pin the software so next year matches" | 2 | workflow-management (digest/lockfile) |
| "Prove a tool change didn't alter my output" | 3 | reproducibility/workflow-testing |
| "Trace every output to its real inputs" | 5 | reproducibility/workflow-provenance |
| "Publish so others can cite and re-run it" | 1+2+3+5 | workflow-provenance (WorkflowHub/Dockstore) |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| "Reproducible" pipeline gives a different result next year | availability mistaken for reproducibility; tools unpinned | pin env by digest/lockfile (pillar 2); this is not optional |
| Reviewer cannot follow the analysis | monolithic script, implicit I/O, why/how undocumented | modularize with explicit contracts; document intent |
| Module runs on your laptop only | hard-coded absolute paths, data mixed into code | parameterize inputs; move data to config/samplesheet |
| Diagram contradicts the code | hand-drawn slide, never updated | regenerate from code (BioFlow-Insight, `--rulegraph`, `-with-dag`) |
| Two workflows "use bwa" but results differ | same free-text name, different tool/version | resolve identity via bio.tools/EDAM; pin the digest |
| Method-comparison result is not defensible | old pipeline never re-run on old data before comparing | reproduce the baseline first, then compare |
| Nobody reuses the published modules | never designed for reuse; not curated | design contracts deliberately; adopt/contribute to a curated registry |

## Related Skills

- reproducibility/workflow-testing - Pillar 3: reference data and regression tests that catch silent breakage
- reproducibility/workflow-provenance - Pillar 5: fine-grained provenance and reproducible publishing
- workflow-management/nextflow-pipelines - Pillar 2 in Nextflow: pin containers by digest, DSL2 module conventions
- workflow-management/snakemake-workflows - Pillar 2 in Snakemake: rules, conda lockfiles, `--rulegraph`
- workflow-management/nf-core-pipelines - Adopt a curated, reviewed pipeline before authoring

## References

- Cohen-Boulakia S, Belhajjame K, Collin O, Chopard J, Froidevaux C, et al. 2017. Scientific workflows for computational reproducibility in the life sciences: status, challenges and opportunities. *Future Gener Comput Syst* 75:284-298.
- Marchment G, Brancotte B, Schmit M, Lemoine F, Cohen-Boulakia S. 2024. BioFlow-Insight: facilitating reuse of Nextflow workflows with structure reconstruction and visualization. *NAR Genom Bioinform* 6(3):lqae092.
- Sebe C, Ferret O, Neveol A, Esmailoghli M, Leser U, Cohen-Boulakia S. 2026. Supporting workflow reproducibility by linking bioinformatics tools across papers and executable code (CoPaLink). *arXiv* 2603.08195.
- Cuevas Villarmin C, Naderi N, Cohen-Boulakia S. 2024. Reproducibility in named entity recognition: a case study analysis. *IEEE eScience*.
- Langer BE, Amaral A, Baudement M-O, et al. 2025. Empowering bioinformatics communities with Nextflow and nf-core. *Genome Biol* 26(1):228.
- Sandve GK, Nekrutenko A, Taylor J, Hovig E. 2013. Ten simple rules for reproducible computational research. *PLoS Comput Biol* 9(10):e1003285.
- Grüning B, Chilton J, Köster J, et al. 2018. Practical computational reproducibility in the life sciences. *Cell Syst* 6(6):631-635.
