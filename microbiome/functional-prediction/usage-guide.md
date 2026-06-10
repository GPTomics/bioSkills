# Functional Prediction - Usage Guide

## Overview

PICRUSt2 predicts community functional POTENTIAL from 16S/ITS amplicon ASVs by phylogenetic interpolation: it places each ASV on a tree of ~20,000 reference genomes and reports the gene content of its nearest sequenced relatives as the community's KO/EC/MetaCyc profile. The single fact that governs every decision: this is PREDICTED potential inferred from who-is-there, never measured gene content and never activity or expression. PICRUSt2 never sequences a functional gene from the sample. Accuracy is entirely a function of how well the community is represented in the reference set - credible in the densely-referenced human gut (predicted-vs-shotgun Spearman ~0.8), collapsing toward a restatement of taxonomy in soil, marine, and novel environments. Because the prediction is a deterministic function of the ASV table, a predicted "functional difference" between groups can be the taxonomic difference re-encoded - so predicted-function analysis is hypothesis-generating, and a functional conclusion needs shotgun or metatranscriptomics.

## Prerequisites

```bash
# PICRUSt2 installs as its own conda environment (pulls EPA-ng, gappa, hmmer, castor)
conda create -n picrust2 -c bioconda -c conda-forge picrust2

# Or the QIIME2 plugin path (into a QIIME2 env) - see qiime2-workflow
conda install -c conda-forge -c bioconda q2-picrust2
```

Conceptual prerequisites:
- Inputs come from upstream amplicon processing: a representative ASV FASTA and an ASV abundance table (BIOM or TSV, samples as columns) from amplicon-processing + taxonomy-assignment.
- The reference (tree, alignment, per-genome trait tables) ships with PICRUSt2 - no separate multi-GB download. The release fixes the reference and therefore the accuracy ceiling.
- NSTI is the mandatory quality gate. The result is uninterpretable without reporting it.
- The result is community POTENTIAL, not activity, and not measured gene content.

## Quick Start

Tell your AI agent what you want to do:
- "Predict KO and MetaCyc potential from my 16S ASV table with PICRUSt2"
- "Run PICRUSt2 and report the NSTI distribution and the fraction of reads dropped"
- "Tell me whether PICRUSt2 is appropriate for my soil samples or whether I should use FAPROTAX"

## Example Prompts

### Prediction
> "I have a representative-sequence FASTA and an ASV abundance table from DADA2. Run PICRUSt2 with the recommended maximum-parsimony hidden-state method and produce KO, EC, and MetaCyc pathway tables."

> "Run the q2-picrust2 full pipeline on my QIIME2 feature table and rep-seqs."

### Quality gating
> "Summarize the NSTI distribution from my PICRUSt2 run and report how many ASVs and what fraction of reads were dropped at the default NSTI cutoff."

> "My mean NSTI is high - is this prediction trustworthy for marine sediment, and what should I use instead?"

### Method choice
> "I want to know whether this community is nitrifying. Should I use PICRUSt2 or FAPROTAX?"

> "I need measured functional gene content, not a prediction - what should I do?"

### Downstream
> "Run compositionally-aware differential abundance on my predicted MetaCyc pathways across two groups and frame the result correctly."

## What the Agent Will Do

1. Confirm the inputs are a representative ASV FASTA plus an abundance table with samples as columns, and identify the environment (gut vs soil/marine/novel).
2. Run `picrust2_pipeline.py` with `--hsp_method mp` and `--max_nsti 2`, producing KO/EC/MetaCyc tables and the NSTI file.
3. Read `marker_predicted_and_nsti.tsv.gz` and report mean/median NSTI, the distribution, and the number of ASVs and fraction of reads dropped at the NSTI cutoff.
4. Flag low reference coverage (high NSTI) and recommend FAPROTAX or shotgun where appropriate.
5. Attach human-readable descriptions to the pathway table with `add_descriptions.py`.
6. For between-group comparisons, hand off to compositionally-aware DA (>=2 tools, report the intersection) and frame predicted differences as hypothesis-generating.
7. Restrict all interpretation to "potential / predicted capacity" - never activity or expression.

## Tips

- Report NSTI every time: mean, median, distribution, and the fraction of reads dropped at `--max_nsti 2`. A run that loses a large read fraction predicted function for a different community than was sampled.
- Treat predicted function as taxonomy re-encoded. Reporting "differed taxonomically AND functionally" double-counts one finding - the two are not independent.
- Restrict claims to potential. "Increased butyrate-production capacity" is defensible; "increased butyrate production" is not.
- PICRUSt2 is strongest in the human gut and weak in soil/marine/novel systems. For biogeochemical-cycle questions in environmental samples, prefer FAPROTAX; for real measured function, use shotgun HUMAnN.
- Use `--hsp_method mp` (maximum parsimony) - the recommended default. `pic` is faster but not recommended.
- The NSTI output file is `marker_predicted_and_nsti.tsv.gz` and the column is `metadata_NSTI`.
- For differential abundance, pass count-like predicted abundances (features as rows) to ALDEx2, not relab-normalized output, and run >=2 compositional tools.
- 16S cannot see strain-level function (accessory genome, HGT, plasmids, pathogenicity islands); the species core is the resolution ceiling regardless of NSTI.

## Related Skills

- amplicon-processing - Generate the ASV table and representative sequences consumed here
- taxonomy-assignment - Taxonomic labels for the same ASVs
- differential-abundance - Compositional DA of the predicted KO/pathway table
- qiime2-workflow - The q2-picrust2 plugin path inside QIIME2
- metagenomics/functional-profiling - MEASURED shotgun function (HUMAnN); the predicted-vs-measured wall
- pathway-analysis/go-enrichment - Reading/enriching the predicted KO/MetaCyc lists
- workflows/microbiome-pipeline - End-to-end amplicon pipeline
