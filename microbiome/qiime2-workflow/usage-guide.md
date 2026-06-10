# QIIME2 Workflow - Usage Guide

## Overview

This is the GLUE skill for the microbiome category. QIIME2's value is not a different algorithm - it wraps DADA2, naive-Bayes/VSEARCH classifiers, MAFFT/FastTree, and the same diversity metrics the R workflow uses. Its value is that a `.qza` artifact is data PLUS its entire executable history: every action, parameter, version, and citation is recorded as provenance, and the semantic-type system refuses to hand the wrong data to an action before it runs. That provenance is the deliverable. The cost is that there is no `cat`-ing the data - one works through the framework or exports and loses the history.

This skill owns the machinery: the artifact/type/provenance model, import (Casava/manifest/EMP/BIOM), export, the Metadata object, provenance replay, and the q2cli vs Artifact API interfaces. It DEFERS every scientific choice to the five sibling skills - denoising parameters to amplicon-processing, classifier/database choice to taxonomy-assignment, sampling depth and diversity metrics to diversity-analysis, differential-abundance tool choice to differential-abundance, and PICRUSt2 to functional-prediction. Shotgun metagenomics uses a different distribution (moshpit) and lives in the metagenomics category.

The release model is a moving target: QIIME2 is calendar-versioned (`YYYY.RELEASE`), ships as separate distributions (amplicon/moshpit/pathogenome/tiny), the framework was rebranded `rachis` in 2026.1, and the `amplicon` distribution is being renamed `qiime2` in 2026.4. Verify the current release and distribution names against the live install page before pinning anything.

## Prerequisites

Install QIIME2 as its own conda environment from the current distribution YAML. The exact filename/URL changes every release - pull it from the live install page (library.qiime2.org / amplicon-docs.qiime2.org), do not hard-code an old URL.

```bash
# 2026.1-era amplicon distribution (verify the current filename on the install page)
conda env create \
  --name qiime2-amplicon-2026.1 \
  --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2026.1-py310-linux-conda.yml
conda activate qiime2-amplicon-2026.1
qiime info     # confirm the release, distribution, and installed plugins
```

Conceptual prerequisites:
- Reads must be demultiplexed FASTQ (or EMP-multiplexed with a barcode column, then demultiplexed in QIIME2).
- Modern Illumina is Phred33 - the offset is baked into the import format name.
- A taxonomy classifier `.qza` must match the QIIME2 release (it is pinned to its scikit-learn version); download the release-namespaced one.
- A metadata TSV with a recognized ID column drives nearly every action; validate it with Keemei before running.
- Reference databases (SILVA/GTDB/UNITE) and pre-trained classifiers are large downloads.

## Quick Start

Tell your AI agent what you want to do:
- "Import my paired-end 16S reads into QIIME2 with a manifest"
- "Check what semantic type this .qza is"
- "Replay the provenance of this artifact to recover the commands"
- "Export my feature table to phyloseq without losing the provenance until the last step"
- "Fix this QIIME2 semantic-type / Phred / classifier-version error"

## Example Prompts

### Import and inspect
> "I have demultiplexed paired-end FASTQ at arbitrary paths. Write a V2 manifest with absolute paths, import them as a SampleData[PairedEndSequencesWithQuality] artifact with the right Phred offset, and summarize the demux to confirm the qualities decoded."

> "What is inside this .qza - its UUID, semantic type, and format? And is the archive valid?"

### Provenance
> "Here is a single .qza a collaborator sent. Replay its provenance to regenerate the commands and pull the citations BibTeX."

> "I want to share my results with collaborators who do not have QIIME2 installed - how do they view the visualizations and the provenance?"

### Orchestration and export
> "Orchestrate the amplicon pipeline through QIIME2 - import, denoise, classify, build a tree, run core diversity, and run ANCOM-BC - but defer the parameter choices to the right skills."

> "Export my feature table, taxonomy, and tree into a phyloseq object for custom R analysis, and tell me where the provenance chain ends."

### Metadata and errors
> "My metadata has an integer subject-ID column. Make sure QIIME2 treats it as categorical, not a continuous covariate, and validate the sheet."

> "My classifier.qza is from a 2024 release and I am on 2026.1 - it throws a scikit-learn version error. What is the fix?"

## What the Agent Will Do

1. Confirm the installed QIIME2 release and distribution (`qiime info`) and verify distribution/plugin names against the current docs.
2. Import raw data as a typed artifact - choosing the on-ramp (Casava directory, V2 manifest, EMP + demux, or BIOM) and the correct Phred offset.
3. Summarize the demux (`qiime demux summarize`) so the truncation decision (owned by amplicon-processing) can be made.
4. Orchestrate the pipeline - denoise, classify, build phylogeny, run core diversity, run ANCOM-BC - calling each action while deferring the scientific parameter choice to the owning sibling skill.
5. Use `qiime tools peek`/`validate` to confirm semantic types before wiring actions together.
6. Annotate the Metadata `#q2:types` for ID/batch columns and recommend Keemei validation.
7. Replay provenance to recover commands and citations from an artifact when reproducing a shared analysis.
8. Export to BIOM/TSV or read into phyloseq only at the last step, flagging the provenance-loss point.

## Tips

- Keep upstream `.qza` files, not just the `.qzv` visualizations - a `.qzv` is terminal and cannot be fed into another action.
- Export is a one-way door: it drops the provenance. Export at the last possible step, or use `qiime2R::qza_to_phyloseq` so the chain survives as far as possible.
- A semantic-type error is the framework working, not a bug. Run `qiime tools peek` and fix the upstream action - do not re-import to coerce the type.
- The classifier is part of the method: match it to the release or retrain it. The `data.qiime2.org/<release>/common/...` URLs are release-namespaced.
- view.qiime2.org renders both visualizations and the provenance DAG client-side with no install - the simplest way to share auditable results.
- The Artifact API (`from qiime2 import Artifact, Metadata`) is better for notebooks and embedding QIIME2 in a larger Python pipeline; provenance is identical to q2cli.
- The legacy q2studio desktop GUI is dead - use Galaxy or view.qiime2.org for no-CLI workflows.
- When the framework overhead is not worth it (a one-off custom R analysis), use qiime2R/export and own the provenance loss honestly rather than fighting the framework.

## Related Skills

- amplicon-processing - DADA2 denoising parameters this skill defers (trunc/trim/maxEE, DADA2 vs Deblur)
- taxonomy-assignment - Classifier and reference-database choice and training behind classify-sklearn
- diversity-analysis - Sampling depth, diversity metric, rarefaction, and the PERMANOVA-vs-dispersion confound
- differential-abundance - DA tool choice and consensus behind composition ancombc
- functional-prediction - PICRUSt2 functional prediction from the feature table
- metagenomics/kraken-classification - Shotgun (moshpit distribution) read classification, not amplicon
- phylogenetics/tree-io - Phylogenetic tree I/O for UniFrac / Faith PD
- read-qc/adapter-trimming - cutadapt primer removal before import/denoising
- workflows/microbiome-pipeline - End-to-end amplicon pipeline
