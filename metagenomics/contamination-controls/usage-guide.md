# Contamination Controls - Usage Guide

## Overview
This skill cleans a shotgun metagenome of everything that is not the target community before profiling: host-read depletion, reagent/kitome contamination control with blanks and decontam, mock-community validation, and depth-adequacy checks. The unifying idea is that a metagenomic result is a position in a choice-chain (extraction, depletion, depth, classifier, database, normalization), not a direct observation - so absence means not-detectable-by-this-chain, and in low biomass the entire community can be reagent contamination. Controls and a consistent chain are what make a result defensible.

## Prerequisites
```bash
conda install -c bioconda hostile nonpareil bowtie2
Rscript -e "BiocManager::install('decontam')"
# Host index (T2T-CHM13-based) downloads on first Hostile run.
```

Conceptual prerequisites:
- Extraction is the experiment: bead-beat, hold one method constant, validate lysis with a whole-cell mock.
- Low-biomass samples need extraction blanks, DNA quantification, and decontam - non-negotiable.
- Remove host reads first (analytical and data-sharing/ethics reasons); report reads removed.
- A confident classifier call can still be wrong if the reference is contaminated.

## Quick Start
Tell your AI agent what you want to do:
- "Remove human reads from my stool metagenome before profiling"
- "Use my extraction blanks to flag reagent contaminants with decontam"
- "Check whether my low-biomass skin samples are real or kitome"
- "Tell me if I sequenced deeply enough to call rare taxa"

## Example Prompts

### Host depletion
> "These are human tissue shotgun reads, mostly host. Deplete human reads with a T2T-CHM13 index, report how many were removed, and check I still have enough microbial depth."

### Kitome control
> "I have 8 low-biomass BAL samples and 2 extraction blanks with DNA concentrations. Run decontam with the prevalence and frequency methods, use the aggressive low-biomass threshold, and show which flagged taxa are canonical kitome genera."

### Reality check
> "My low-biomass sample shows a novel Bradyrhizobium-dominated community. Help me decide whether that is real or reagent contamination."

### Depth adequacy
> "Run Nonpareil and tell me whether my depth supports a claim that a pathogen is absent."

## What the Agent Will Do
1. Host-deplete with a T2T-CHM13 index and a high-sensitivity aligner, reporting reads removed.
2. Run decontam on the classifier output table using blanks and DNA concentration, per batch.
3. Inspect flagged contaminants against the canonical kitome genera before removing them.
4. Validate the pipeline with a mock community when limit of detection or lysis is in question.
5. Check depth with Nonpareil and refuse to interpret absence below the limit of detection.
6. Keep the whole chain consistent and report every link (kit/lot, host index, blanks, mock version, depth).

## Tips
- The lower the biomass, the larger the kitome fraction - at near-zero biomass the signal is the kitome.
- Use a whole-cell mock to test extraction/lysis and a DNA mock to test the classifier/library.
- Prefer T2T-CHM13 over GRCh38 and mask host rDNA so you do not delete real microbial reads.
- decontam runs on the feature table (Bracken/MetaPhlAn output), not raw reads; run it per batch.
- Do not meta-analyze across studies that differ in extraction, depth, or database.

## Control Types

| Control | Tests | Where |
|---------|-------|-------|
| Extraction/no-template blank | the kitome for this batch/lot | full workflow, >=1 per batch |
| Whole-cell mock (ZymoBIOMICS) | extraction/lysis bias + classifier | alongside samples |
| DNA mock | classifier/library accuracy + limit of detection | alongside samples |
| DNA concentration (Qubit/PicoGreen) | the decontam frequency signal | per sample |

## Resources
- [decontam](https://benjjneb.github.io/decontam/)
- [Hostile](https://github.com/bede/hostile)
- [Nonpareil](https://nonpareil.readthedocs.io/)

## Related Skills

- kraken-classification - Classification after host removal; database bias
- metaphlan-profiling - Marker-gene profiling after cleanup
- abundance-estimation - decontam runs on the abundance table
- metagenome-visualization - Plot blanks alongside samples
- read-qc/adapter-trimming - Generic trimming before this step
- genome-assembly/metagenome-assembly - MAG-level decontamination
- workflows/metagenomics-pipeline - End-to-end pipeline with a controls stage
