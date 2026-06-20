# Contamination Screening - Usage Guide

## Overview
Contamination comes in five distinct classes, each with its own detection method and fix. A species screen (FastQ Screen, Kraken2) maps reads to genome references and answers "what organisms are here?" -- but it is structurally blind to same-species cross-sample contamination and sample swaps, which need SNP fingerprinting (verifyBamID2, NGSCheckMate, somalier, conpair). Any human or single-species cohort needs both. The default action is to report contamination, not filter, because removing reads biases composition.

## Prerequisites
```bash
conda install -c bioconda fastq-screen kraken2 bracken bbmap sortmerna somalier
# Pre-built FastQ Screen databases
fastq_screen --get_genomes
```

## Quick Start
Tell your AI agent what you want to do:
- "Screen my FASTQ files for cross-species contamination"
- "Check if my human sample has mouse (PDX host) contamination"
- "Verify my human samples are not swapped or mixed"
- "Estimate the contamination fraction in my BAM"

## Example Prompts

### Cross-species screening
> "Run FastQ Screen on my samples against human, mouse, and common contaminants"

> "Classify my reads with Kraken2 to find a bacterial contaminant"

### Same-species (swaps and mixtures)
> "My RNA-seq results look wrong; check whether samples are swapped"

> "Estimate within-species contamination in my tumor and normal BAMs"

### PDX and remediation
> "Separate human graft from mouse host reads in my PDX sample"

> "Remove PhiX spike-in from my reads before assembly"

## What the Agent Will Do
1. Identify which contamination class is in question (cross-species, swap, vector, rRNA, cell-line)
2. Run the matching detector (genome panel for organisms; SNP fingerprints for individuals)
3. Read the FastQ Screen bar chart by category (unexpected one-hit-one-genome = contamination; multiple-genomes = homology)
4. Recommend report vs filter, preferring combined-reference alignment for PDX
5. For human cohorts, add a SNP-fingerprint identity/contamination check

## Interpretation Guide

### FastQ Screen bar chart (cross-species only)
- An unexpected genome with high One_hit_one_genome = real foreign DNA.
- Reads spread across multiple genomes' multiple-hits categories = conserved/homologous sequence (rRNA, mito), not contamination.
- High Hit_no_genomes = adapter dimer, a missing reference, or a novel organism.

### Same-species signals (NOT visible to FastQ Screen)
- verifyBamID2 FREEMIX > ~0.02 = within-species contamination.
- FREEMIX near 0 but CHIPMIX near 1 = a sample swap, not contamination.
- somalier/NGSCheckMate off-diagonal identity = swapped or mixed samples.

## Tips
- Run an organism screen early, before investing in alignment.
- A species screen cannot detect a human-A/human-B swap; use SNP fingerprints for that.
- Index hopping on patterned flowcells creates phantom low-VAF variants; use unique dual indexing for ctDNA/single-cell/somatic.
- Default to report-only; filter only a NAMED contaminant (PhiX, adapters) with a precise k-mer remover (BBDuk).
- For PDX, prefer aligning to a combined human+mouse reference over hard pre-filtering.
- Reference-genome contamination can corrupt the screen; trust a curated database and treat surprising single-source hits as artifacts until ruled out.

## Resources
- [FastQ Screen Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
- [Kraken2 / Bracken](https://github.com/DerrickWood/kraken2)
- [somalier](https://github.com/brentp/somalier)

## Related Skills
read-qc/quality-reports - Bimodal GC and overrepresented sequences flag contamination
read-qc/adapter-trimming - Remove adapter contamination
read-qc/rnaseq-qc - rRNA fraction as a prep-efficiency metric
metagenomics/kraken-classification - Deeper taxonomic classification and profiling
variant-calling/joint-calling - Where SNP-fingerprint sample swaps do the most damage
