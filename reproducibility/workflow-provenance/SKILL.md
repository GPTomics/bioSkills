---
name: bio-reproducibility-workflow-provenance
description: Captures the recorded lineage of a bioinformatics workflow run and packages it for reproducible sharing, distinguishing provenance (what produced each output) from reproducibility (the ability to regenerate it). Use when deciding fine-grained why-provenance (link each output only to the inputs that produced it) versus coarse capture that causes provenance overload; modelling lineage with W3C PROV (Entity/Activity/Agent, used/wasGeneratedBy/wasDerivedFrom); choosing a provenance artifact (Workflow Run RO-Crate via nf-prov or runcrate, CWLProv via cwltool --provenance, BioCompute Object for FDA/regulatory submissions, or engine-native Nextflow/Snakemake traces); capturing the environment (container digests, versions) so a re-run can be shown to match; and publishing a runnable, citable workflow (WorkflowHub DOIs, Dockstore, GA4GH WES/TRS/DRS). Routes environment pinning and expected-output testing to the sibling skills.
tool_type: mixed
primary_tool: ro-crate-py
---

## Version Compatibility

Reference examples tested with: ro-crate-py (rocrate) 0.10+, runcrate 0.5+, cwltool 3.1+, Nextflow 24.04+, nf-prov 1.2+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show rocrate runcrate cwltool` then `python -c "help('rocrate')"`
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and
adapt to the actual API rather than retrying.

Note: the RO-Crate profiles evolve (Workflow RO-Crate for the workflow, Provenance/Workflow Run
RO-Crate for a specific execution); confirm the current profile URIs. `cwltool --provenance <dir>`
emits a CWLProv Research Object; nf-prov emits Workflow Run RO-Crate and BioCompute Object formats.

# Workflow Provenance

**"Record where every output came from, and package the run so others can re-run and cite it"** -> Emit machine-readable lineage (which activity, using which inputs and which tool version, produced which output) and bundle the workflow, its metadata, and that lineage into a portable, DOI-citable artifact. Provenance is the record; reproducibility is the ability to act on it.

## The governing principle: provenance != reproducibility, and fine-grained beats overload

Provenance is the recorded lineage of every output; reproducibility is the ability to regenerate it (workflow-management README; Cohen-Boulakia 2017 *Future Gener Comput Syst* 75:284-298). You can have detailed provenance of a run you cannot reproduce (the environment was not pinned), and you can reproduce a run whose provenance you never recorded. Both are needed, and they are captured by different means.

The research-grade hard part is granularity. Naively, "re-run and track everything at each step" records that every output depends on every input of its step. But most outputs depend only on *some* inputs (a step that takes a list and returns the difference of its last two elements did not use the rest). Recording all-to-all produces **provenance overload**: so much lineage that neither a tool nor a biologist can see what actually caused what. The valuable goal is **fine-grained (why-)provenance**: each output linked only to the inputs that genuinely produced it. The analogy is a social network studied while keeping different link types distinct (friendships vs phone calls); collapsing them destroys the readable signal.

| Property | Provenance | Reproducibility |
|----------|-----------|-----------------|
| Question | what produced this output? | can I regenerate this output? |
| Artifact | PROV graph, RO-Crate, BCO, trace | pinned env + reference data + tests |
| Failure mode | overload (too coarse) or gaps (uncaptured) | drift (unpinned) |
| Owned by | this skill | workflow-management (env) + workflow-testing |

## The W3C PROV model (the vocabulary underneath every format)

PROV-DM has three node types and a handful of edges; RO-Crate, CWLProv, and BCO all serialize this.

| PROV term | Bioinformatics meaning | Key relations |
|-----------|------------------------|---------------|
| Entity | a data object: FASTQ, BAM, VCF, a parameter file | `wasGeneratedBy`, `wasDerivedFrom` (entity->entity) |
| Activity | a step execution: an alignment run | `used` (activity->entity), `wasInformedBy` |
| Agent | a person, an organization, or the software/container | `wasAssociatedWith` (activity->agent), `wasAttributedTo` |

Fine-grained provenance is exactly the discipline of drawing `used` and `wasDerivedFrom` edges only where a real dependency exists, rather than a complete bipartite graph per step.

## Decision: which provenance artifact

| Artifact | Scope | Captures | Emit with | Choose when |
|----------|-------|----------|-----------|-------------|
| Engine-native trace | one run | per-task resource, exit, hash, DAG (not a portable standard) | `nextflow -with-trace -with-report -with-timeline -with-dag`; `snakemake --report` | quick audit/debug of a single run |
| Workflow RO-Crate | the workflow | code + metadata + inputs, portable JSON-LD | ro-crate-py; `nf-core` ships one | packaging a workflow to share/register |
| Workflow Run / Provenance RO-Crate | one execution | the above + PROV of that run | nf-prov plugin; `runcrate` (converts CWLProv) | sharing a specific reproducible run |
| CWLProv Research Object | one CWL run | PROV + inputs/outputs + intermediate data | `cwltool --provenance ro/ ...` | CWL pipelines; strongest out-of-box provenance |
| BioCompute Object (IEEE 2791) | a pipeline for regulators | domains (usability/description/execution/parametric/io/error) | nf-prov (bco format); BCO tools | FDA/clinical submissions, regulated contexts |

CWL ships the strongest standardized provenance artifact by default; RO-Crate is the portable packaging layer the community is converging on (Soiland-Reyes 2022 *Data Science* 5:97-138); Workflow Run RO-Crate is the profile for a specific execution.

## Engine-native trace (the cheap first layer)

```groovy
// nextflow.config -- turn on run-level provenance/trace and a portable Run RO-Crate.
trace   { enabled = true; file = 'prov/trace.txt';   fields = 'task_id,hash,name,status,container,realtime' }
report  { enabled = true; file = 'prov/report.html' }         // resources + versions per task
dag     { enabled = true; file = 'prov/dag.mmd' }             // structure, code-derived
timeline{ enabled = true; file = 'prov/timeline.html' }
plugins { id 'nf-prov' }                                       // emits Workflow Run RO-Crate / BCO
prov    { formats { wrroc { file = 'prov/ro-crate-metadata.json' } } }
```

The trace records the container reference and task hash per step, which is the link between provenance (what ran) and reproducibility (pin that container by digest so it can run again). A `:latest` tag in the trace is a provenance record of an *unreproducible* run.

## Build and read an RO-Crate (portable packaging)

```python
# Reference: rocrate 0.10+ | Verify API if version differs
# Package a workflow + its metadata into a Workflow RO-Crate that others can register/run.
from rocrate.rocrate import ROCrate
from rocrate.model.person import Person

crate = ROCrate()
wf = crate.add_workflow("main.nf", main=True, lang="nextflow",
                        properties={"name": "Variant calling", "version": "1.2.0"})
crate.add(Person(crate, "#author", {"name": "A. Researcher", "affiliation": "Lab"}))
crate.metadata.write("crate/")          # writes ro-crate-metadata.json (JSON-LD)
```

```bash
# Convert a CWLProv Research Object into a Provenance Run Crate (records one execution).
cwltool --provenance ro/ workflow.cwl inputs.yml     # produce the CWLProv RO
runcrate convert ro/ --output run-crate/              # -> Workflow Run RO-Crate
runcrate report run-crate/                            # human-readable lineage of the run
```

## Decision: publishing a runnable, citable workflow

| Goal | Platform / API | What it gives |
|------|----------------|---------------|
| Register + get a DOI to cite | WorkflowHub (RO-Crate based) | versioned, citable workflow entry (Gustafsson 2024) |
| Discoverable, runnable across platforms | Dockstore (GA4GH TRS) | CWL/WDL/Nextflow/Galaxy, TRS registration |
| Search before authoring | Dockstore / WorkflowHub | adopt an existing curated workflow (see by-design) |
| Run the same workflow on any compliant platform | GA4GH WES | portable execution API |
| Reference tools / data by standard ID | GA4GH TRS / DRS | resolvable tool and data identifiers |

Search these registries before authoring: a curated, reviewed workflow is reused far more than a bespoke one and already carries provenance and tests.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Provenance graph is unreadable, everything depends on everything | coarse all-inputs-to-all-outputs capture | record fine-grained `used`/`wasDerivedFrom` only for real dependencies |
| Detailed provenance but the run will not reproduce | environment never pinned; trace shows `:latest` | pin container by digest / conda lockfile (workflow-management) |
| RO-Crate has the workflow but no run lineage | used the Workflow RO-Crate profile, not the Run profile | emit a Workflow Run / Provenance Run Crate (nf-prov, runcrate) |
| Lineage lost after publishing outputs | `publishDir` symlinks broke the derivation trail | capture provenance from `work/`/the RO, not from published symlinks |
| Trace full of PID/timestamp/host noise, no stable identity | recorded volatile fields as identity | key on content/task hash and container digest, not PID/mtime |
| Cannot cite the workflow in a paper | published to a repo with no DOI | register on WorkflowHub for a versioned DOI |
| Regulatory reviewer cannot follow the pipeline | free-text methods, no structured provenance | emit a BioCompute Object (IEEE 2791) |
| Re-run diverges and no one knows which step changed | no per-task hash/version record | keep the engine trace (hash + container) with the run |

## Related Skills

- reproducibility/by-design - Provenance is pillar 5; provenance vs reproducibility framing
- reproducibility/workflow-testing - Records which tool versions produced a tested/expected output
- workflow-management/cwl-workflows - CWLProv (`--provenance`), the strongest out-of-box provenance artifact
- workflow-management/nextflow-pipelines - `-with-trace`/`-with-report`, nf-prov, and pinning containers by digest
- workflow-management/nf-core-pipelines - Curated pipelines that ship an RO-Crate and version metadata

## References

- Cohen-Boulakia S, Belhajjame K, Collin O, et al. 2017. Scientific workflows for computational reproducibility in the life sciences: status, challenges and opportunities. *Future Gener Comput Syst* 75:284-298.
- Soiland-Reyes S, Sefton P, Crosas M, et al. 2022. Packaging research artefacts with RO-Crate. *Data Sci* 5(2):97-138.
- Khan FZ, Soiland-Reyes S, Sinnott RO, Lonie A, Goble C, Crusoe MR. 2019. Sharing interoperable workflow provenance: a review of best practices and their implementation in CWLProv. *GigaScience* 8(11):giz095.
- Goble C, Cohen-Boulakia S, Soiland-Reyes S, et al. 2020. FAIR computational workflows. *Data Intelligence* 2(1-2):108-121.
- Gustafsson OJR, Wilkinson SR, Bacall F, et al. 2024. WorkflowHub: a registry for computational workflows. *arXiv* 2410.06941.
- Alterovitz G, Dean D, Goble C, et al. 2018. Enabling precision medicine via standard communication of HTS provenance, analysis, and results (BioCompute Object). *PLoS Biol* 16(12):e3000099.
