# Workflow Testing - Usage Guide

## Overview

Workflow testing ships a sample input and its expected output alongside a bioinformatics pipeline so anyone can verify it still works in minutes, and so continuous integration flags the moment a tool or container version change silently alters a result. Providing reference data and reference results is the community's weakest reproducibility practice (Cohen-Boulakia 2017), and a pipeline that merely exits 0 has proved nothing about correctness. This skill covers the three test tiers, generating shareable synthetic data without disclosing sensitive samples, and assertion strategies for outputs that are not bytewise reproducible.

## Prerequisites

```bash
# Test frameworks (choose per engine)
curl -fsSL https://get.nf-test.com | bash      # nf-test (Nextflow; needs a JVM)
pip install pytest-workflow                     # engine-agnostic YAML tests
pip install cwltest                             # CWL

# Test-data simulation and downsampling
conda install -c bioconda dwgsim art insilicoseq seqtk samtools
# RNA-seq truth data:  R -e 'BiocManager::install("polyester")'
```

Inputs: a runnable pipeline, a tiny reference dataset (curated, downsampled, or simulated with a fixed seed), and a known-good expected output to snapshot.

## Quick Start

Tell your AI agent what you want to do:
- "Add a stub/dry-run wiring test and a tiny smoke test to my pipeline"
- "Generate synthetic reads with known variants so I can test without patient data"
- "Write assertions that normalize the VCF before comparing so header timestamps do not fail the test"
- "Set up nf-test snapshot testing and a GitHub Actions matrix over docker and singularity"
- "Make CI alert me when a container version bump changes my output"

## Example Prompts

### Choosing Test Tiers
> "Set up the three tiers for my Nextflow pipeline: a stub wiring test, a tiny smoke test that runs the tools, and a scheduled full-dataset run."

> "My CI is too slow; move the full dataset out of the per-PR tests and keep the smoke data tiny."

### Generating Test Data
> "Simulate 500 paired-end reads from chr20 with two known SNVs using dwgsim and a fixed seed, so I can test variant calling without real samples."

> "Downsample this public 1000G BAM to 1% with a recorded seed to make a shareable fixture."

### Assertion Strategy
> "My md5 test fails every run even though the VCF looks identical; fix it to compare records without the header."

> "Add snapshot testing to my nf-test so a silent output change shows up as a diff in the PR."

### CI and Alerting
> "Write a GitHub Actions workflow that runs stub then nf-test over a docker/singularity matrix on every PR."

> "Pin my containers by digest and make the smoke test go red when someone bumps a version."

## What the Agent Will Do

1. Pick the test tiers (stub wiring, tiny smoke, scheduled integration) and gate the expensive ones behind the cheap ones.
2. Generate or downsample shareable reference data with a fixed seed, keeping it CI-sized and versioned with the tests.
3. Choose an assertion appropriate to the output type instead of a naive whole-file md5 that fails on timestamps and ordering.
4. Set up snapshot testing so an intentional behavior change is a reviewable diff and a silent one turns the build red.
5. Select the framework (nf-test, pytest-workflow, snakemake unit tests, cwltest, miniwdl) matching the engine.
6. Wire CI that runs the cheap tiers on every pull request and the full run on a schedule.
7. Pin containers by digest so a version bump is caught by the snapshot diff rather than discovered mid-experiment.

## Tips

- A pipeline exiting 0 proves only that nothing crashed; add content and structural assertions, not just an exit-code check.
- Genomic outputs are rarely bytewise reproducible; strip or normalize headers and sort records before comparing VCF/BAM/tables.
- Always fix the simulator RNG seed and record it in the fixture, or the "reference" data is not itself reproducible.
- Keep smoke-tier data tiny enough to run on every PR; put full reference datasets in a scheduled job.
- Treat a snapshot change in a pull request as a claim that behavior changed on purpose, and review the diff.
- Mirror the nf-core pattern: tiny inputs wired through `-profile test`, so users can smoke-test with one command.
- A container/tool version bump that changes results while every step still exits 0 is exactly the failure a regression snapshot catches; pin by digest so the diff appears.

## Related Skills

- reproducibility/by-design - Why testing is pillar 3 of reproducibility by design
- reproducibility/workflow-provenance - Records which tool versions produced a tested output
- workflow-management/nf-core-pipelines - The `-profile test` and `nf-core/test-datasets` conventions
- workflow-management/nextflow-pipelines - Stub wiring tests and pinning containers by digest
- read-qc/quality-reports - QC outputs that make good structural assertion targets
