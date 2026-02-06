# ncRNA Search - Usage Guide

## Overview
Search for non-coding RNA homologs and classify RNA families using Infernal covariance model (CM) searches against the Rfam database. Covariance models capture both sequence and secondary structure conservation, providing higher sensitivity for structured RNAs than BLAST or profile HMMs alone.

## Prerequisites
```bash
# Infernal
conda install -c bioconda infernal

# Download and prepare Rfam database
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
cmpress Rfam.cm

# Clan information for overlap resolution
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin

# Python dependencies
pip install biopython pandas
```

## Quick Start
Tell your AI agent what you want to do:
- "Search my RNA sequence against Rfam to identify what family it belongs to"
- "Scan my transcriptome assembly for known ncRNA families"
- "Build a covariance model for my novel RNA family"
- "Find all tRNAs in my genome assembly"
- "Classify non-coding transcripts from my RNA-seq experiment"

## Example Prompts

### Rfam Classification
> "I have a set of non-coding transcripts. Scan them against Rfam to classify by family."

> "Search this single RNA sequence against Rfam and tell me what it is."

### Genome-wide ncRNA Search
> "Find all rRNAs, tRNAs, and snoRNAs in my bacterial genome using Infernal."

> "Search for a specific Rfam family (RF00005, tRNA) in my genome assembly."

### Custom Covariance Models
> "I have a Stockholm alignment of a novel RNA family with consensus structure. Build a CM and search my genome."

> "Take the hits from my initial search and refine the covariance model iteratively."

### Parsing and Analysis
> "Parse the Infernal tabular output and summarize ncRNA family assignments."

> "Extract the sequences of all significant Rfam hits from my FASTA file."

## What the Agent Will Do
1. Verify Rfam database is downloaded and indexed (cmpress)
2. Run cmscan (query vs Rfam) or cmsearch (CM vs database) with appropriate thresholds
3. Apply clan overlap resolution to remove redundant family hits
4. Parse tabular output and filter by E-value or gathering threshold
5. Summarize results by family, count, and genomic location

## Tips
- **Gathering thresholds** - Use `--cut_ga` for Rfam searches; these are curated per-family cutoffs that balance sensitivity and specificity
- **Clan overlap** - Always use `--clanin` and `--oclan` with Rfam to resolve overlapping hits from related families
- **Custom CMs** - Always run `cmcalibrate` after `cmbuild`; without calibration, E-values are unreliable
- **Speed** - Use `--noali` to skip alignment output for large-scale searches; add `--rfam` for accelerated heuristics
- **Truncated hits** - Hits at contig/sequence boundaries may be incomplete; check the trunc column in output

## Related Skills
- secondary-structure-prediction - Predict structures for novel ncRNA candidates
- genome-annotation/ncrna-annotation - Genome-wide ncRNA annotation pipelines
- alignment/msa-statistics - Evaluate alignment quality before CM building
