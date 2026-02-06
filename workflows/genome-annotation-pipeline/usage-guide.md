# Genome Annotation Pipeline - Usage Guide

## Overview
Complete genome annotation workflow from assembled contigs to functional annotation. Supports both prokaryotic (Bakta one-step) and eukaryotic (RepeatMasker, BRAKER3, eggNOG-mapper, Infernal) paths with QC checkpoints at each stage.

## Prerequisites
```bash
# Prokaryotic path
conda install -c bioconda bakta
bakta_db download --output /path/to/bakta_db --type full

# Eukaryotic path
conda install -c bioconda repeatmodeler repeatmasker braker3 busco agat
conda install -c bioconda eggnog-mapper interproscan trnascan-se infernal

# Python dependencies
pip install pandas biopython
```

**Input data:**
- Assembled genome FASTA (contigs or scaffolds)
- For eukaryotes: RNA-seq BAM and/or protein evidence FASTA (e.g., OrthoDB)

## Quick Start
Tell your AI agent what you want to do:
- "Annotate my newly assembled bacterial genome"
- "Run a full eukaryotic genome annotation pipeline"
- "Annotate my fungal assembly with repeat masking and gene prediction"
- "Set up a prokaryotic annotation with Bakta"

## Example Prompts

### Prokaryotic
> "I have an assembled E. coli genome. Run Bakta annotation with NCBI-compatible locus tags."

> "Annotate my bacterial contigs and check if the gene count looks reasonable."

### Eukaryotic
> "I have an assembled insect genome with RNA-seq alignments. Run the full annotation pipeline: repeats, gene prediction, functional annotation."

> "Mask repeats in my plant genome, predict genes with BRAKER3, and annotate them with eggNOG-mapper."

### QC Focus
> "Check whether my genome annotation has reasonable BUSCO completeness and gene counts."

> "My annotation only has 40% functional assignment. Help me troubleshoot."

## What the Agent Will Do
1. Assess genome type (prokaryotic or eukaryotic)
2. For prokaryotes: run Bakta one-step annotation
3. For eukaryotes: mask repeats with RepeatModeler/RepeatMasker
4. Predict genes with BRAKER3 using available evidence
5. Run BUSCO to check annotation completeness
6. Assign functions with eggNOG-mapper and InterProScan
7. Annotate ncRNAs with tRNAscan-SE and Infernal
8. Merge annotations and validate final GFF3
9. Report QC metrics at each stage

## Tips
- Always soft-mask (`-xsmall`) rather than hard-mask for gene prediction input
- BRAKER3 with both RNA-seq and protein evidence gives best results
- Use the full Bakta database for comprehensive prokaryotic annotation
- BUSCO completeness > 90% is the target for well-assembled genomes
- If functional annotation is < 60%, try broader taxonomy scope in eggNOG
- Check repeat content against expectations for the taxon before proceeding
- Run QUAST and BUSCO on the assembly itself before starting annotation

## Related Skills
- genome-annotation/prokaryotic-annotation - Bakta and Prokka details
- genome-annotation/eukaryotic-gene-prediction - BRAKER3 and AUGUSTUS options
- genome-annotation/functional-annotation - eggNOG-mapper and InterProScan
- genome-assembly/assembly-qc - Pre-annotation assembly quality checks
- genome-intervals/gtf-gff-handling - GFF3/GTF parsing and manipulation
