# Workflow Provenance - Usage Guide

## Overview

Workflow provenance is the recorded lineage of a run: which activity, using which inputs and which tool version, produced which output. It is distinct from reproducibility (the ability to regenerate the output) and both are needed. This skill covers fine-grained why-provenance (linking each output only to the inputs that genuinely produced it, avoiding provenance overload), the W3C PROV model that underlies every format, the choice of provenance artifact (RO-Crate, CWLProv, BioCompute Object, engine-native traces), and publishing a runnable, citable workflow on WorkflowHub or Dockstore.

## Prerequisites

```bash
# RO-Crate build, Provenance Run Crate conversion, and CWLProv
pip install rocrate runcrate cwltool

# Engine-native provenance
curl -s https://get.nextflow.io | bash     # Nextflow -with-trace/-with-report; nf-prov plugin
pip install snakemake                        # Snakemake --report
```

Inputs: a completed workflow run (its `work/` directory or a CWLProv Research Object), the workflow code, and container digests/versions for the environment record.

## Quick Start

Tell your AI agent what you want to do:
- "Turn on run-level provenance and a trace for my Nextflow pipeline"
- "Produce a CWLProv Research Object for my CWL run"
- "Package this run as a Workflow Run RO-Crate with fine-grained lineage"
- "Emit a BioCompute Object for a regulatory submission"
- "Publish my workflow to WorkflowHub with a DOI so it can be cited"

## Example Prompts

### Capturing Provenance
> "Configure nextflow.config to emit a trace, an execution report, a DAG, and a Workflow Run RO-Crate via nf-prov."

> "Run my CWL workflow with cwltool --provenance and convert the result to a Provenance Run Crate with runcrate."

### Fine-Grained vs Overload
> "My provenance graph shows every output depending on every input; make it fine-grained so each output links only to the inputs that produced it."

> "Explain the difference between provenance and reproducibility for my pipeline and which artifacts cover each."

### Choosing an Artifact
> "Which provenance artifact should I use for an FDA submission versus for sharing on WorkflowHub?"

> "Build a Workflow RO-Crate that packages my main.nf with author and version metadata."

### Publishing
> "Register my workflow on WorkflowHub for a citable DOI and make it runnable as an RO-Crate."

> "Search Dockstore and WorkflowHub for an existing curated pipeline before I author my own."

## What the Agent Will Do

1. Separate provenance (the lineage record) from reproducibility (the pinned environment plus tests) and cover each with the right artifact.
2. Record fine-grained `used`/`wasDerivedFrom` edges only where a real dependency exists, avoiding provenance overload.
3. Turn on engine-native traces (Nextflow report/trace/DAG, Snakemake report) as the cheap first layer, keyed on task hash and container digest.
4. Choose the artifact by goal: RO-Crate to package/share, CWLProv for CWL runs, BioCompute Object for regulated submissions.
5. Build or convert a Workflow Run RO-Crate (ro-crate-py, runcrate) that bundles code, metadata, and run lineage.
6. Publish to WorkflowHub for a versioned DOI or Dockstore for GA4GH TRS registration, and search these before authoring.
7. Ensure the environment is pinned (routing to workflow-management) so the recorded provenance is of a run that can actually be reproduced.

## Tips

- Provenance and reproducibility are different: you can have detailed lineage of a run you cannot regenerate, and vice versa.
- Fine-grained provenance links each output only to the inputs that produced it; recording all-inputs-to-all-outputs creates unreadable overload.
- RO-Crate, CWLProv, and BioCompute Objects all serialize the same W3C PROV vocabulary (Entity/Activity/Agent).
- The engine trace records the container reference per task; a `:latest` in the trace is a provenance record of an unreproducible run, so pin by digest.
- CWL ships the strongest out-of-box provenance (`cwltool --provenance`); RO-Crate is the portable packaging layer the community is converging on.
- Use the Workflow Run (Provenance) RO-Crate profile for a specific execution, not the plain Workflow RO-Crate profile, when you need run lineage.
- Register on WorkflowHub for a citable DOI, and search WorkflowHub/Dockstore before authoring, since curated workflows are reused far more.

## Related Skills

- reproducibility/by-design - Provenance is pillar 5 of reproducibility by design
- reproducibility/workflow-testing - Records which tool versions produced a tested output
- workflow-management/cwl-workflows - CWLProv, the strongest out-of-box provenance artifact
- workflow-management/nextflow-pipelines - Engine traces, nf-prov, and pinning containers by digest
- workflow-management/nf-core-pipelines - Curated pipelines shipping an RO-Crate and version metadata
