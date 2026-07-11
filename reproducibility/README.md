# reproducibility

## Overview

Reproducibility by design for scientific workflows. The premise, following Cohen-Boulakia (2017), is that putting code on GitHub with a Docker image is *availability*, not reproducibility, and reproducibility is the precondition for reuse: nobody can reuse an analysis they cannot re-run and understand. These skills operationalize three pillars that the workflow-engine skills do not cover on their own: designing modular, comprehensible workflows with explicit tool identity; providing reference data and expected outputs so a workflow can be tested and a silent breakage caught; and capturing fine-grained provenance so every output is traceable to the inputs that actually produced it. Environment pinning (containers by digest, conda by lockfile) is the fourth pillar and lives in workflow-management.

**Tool type:** mixed | **Primary tools:** BioFlow-Insight, nf-test, pytest-workflow, ro-crate-py, runcrate, cwltool

## Skills

| Skill | Description |
|-------|-------------|
| by-design | The doctrine and router: the three levels of reproducibility, availability vs reproducibility vs reuse, the five pillars, modular design with explicit I/O contracts, tool identity/naming (bio.tools/EDAM, CoPaLink), and workflow comprehension/decay (BioFlow-Insight) |
| workflow-testing | Reference and synthetic test data (read simulators, downsampling) without disclosing sensitive data, assertion strategies for non-bytewise outputs (normalize-then-compare, snapshots, tolerance), and CI that alerts when a tool change silently breaks the expected output (nf-test, pytest-workflow, snakemake unit tests, cwltest) |
| workflow-provenance | Fine-grained (why-)provenance vs provenance overload, the W3C PROV model, provenance artifacts (Workflow Run RO-Crate, CWLProv, BioCompute Object, engine-native traces), and reproducible sharing/publishing (WorkflowHub, Dockstore, GA4GH WES/TRS/DRS, DOIs) |

## The three levels of reproducibility

The word "reproducible" hides three different claims (Cohen-Boulakia 2017), and a skill or a reviewer request should always name which one is meant.

| Level | Claim | What it captures | How to test it |
|-------|-------|------------------|----------------|
| L1 computational | Same data + same code + same environment -> same result | Maximal context: why and under what conditions a result was obtained | Re-run in the pinned container on the reference data; byte- or record-compare (workflow-testing) |
| L2 robustness | Vary the method a little; understand what the result depends on | Which parameters/tools matter, which do not | Parameter sweeps, alternative tools on the same data; does the finding move? |
| L3 scientific | Any reasonable method, any small choice -> the same finding holds | The biological conclusion, independent of tooling | Independent re-analysis, meta-analysis across cohorts/pipelines |

L1 is the least interesting on its own but the most demanding to document, and it is the floor the other two stand on. These skills target L1-by-design and give the artifacts (tests, provenance) that make L2 and L3 tractable.

## Example Prompts

- "Restructure my monolithic pipeline into modules with explicit input/output contracts so others can reuse the blocks"
- "Reconstruct and visualize my Nextflow workflow's structure so I can see the steps"
- "Generate tiny synthetic test data and expected outputs so my workflow can be regression-tested in CI"
- "Write nf-test snapshot tests that normalize the VCF before comparing so header timestamps do not cause false failures"
- "Package a completed pipeline run as a Workflow Run RO-Crate with fine-grained provenance"
- "Publish my workflow to WorkflowHub with a DOI and a runnable RO-Crate"

## Requirements

```bash
# Workflow comprehension / structure reconstruction
pip install bioflow-insight

# Testing frameworks (choose per engine)
# nf-test (Nextflow): https://www.nf-test.com  (needs a JVM)
curl -fsSL https://get.nf-test.com | bash
pip install pytest-workflow            # engine-agnostic, YAML-defined
pip install cwltest                    # CWL conformance/regression

# Test-data simulation
conda install -c bioconda dwgsim art insilicoseq seqtk samtools

# Provenance and packaging
pip install rocrate runcrate cwltool   # RO-Crate build, Provenance Run Crate, CWLProv
```

## Related Skills

- **workflow-management** - The engines these practices harden (Snakemake/Nextflow/CWL/WDL/nf-core); the fourth pillar (environment pinning by digest/lockfile) lives here
- **workflows** - End-to-end domain pipelines that should ship with tests and provenance
- **reporting** - Reproducible reports (RMarkdown/Quarto/Jupyter) that consume a workflow's outputs
- **experimental-design** - Design-stage reproducibility (randomization, batch structure, power) upstream of the workflow
