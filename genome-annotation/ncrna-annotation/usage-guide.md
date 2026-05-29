# Non-Coding RNA Annotation - Usage Guide

## Overview

Identify and annotate non-coding RNAs (tRNAs, rRNAs, snoRNAs, snRNAs, riboswitches, regulatory sRNAs) in genome assemblies with Infernal covariance-model search against Rfam, tRNAscan-SE 2.0 for tRNA, barrnap for rRNA, and ARAGORN for tmRNA. The animating principle: ncRNA homology is conserved at the level of base-paired secondary structure, not primary sequence, so covariance models (not BLAST) are the engine - and a homology annotation is a recall floor that is blind to fast-evolving and lineage-specific RNAs. miRNA and lncRNA live in expression/transcript space, not genomic search.

## Prerequisites

```bash
# Infernal + Rfam (record the Rfam release)
conda install -c bioconda infernal
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz && cmpress Rfam.cm
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin

# Specialists
conda install -c bioconda trnascan-se barrnap aragorn

# Python utilities
pip install pandas biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "Find all conserved-family non-coding RNAs in my genome assembly"
- "Annotate tRNAs in my bacterial genome and report the high-confidence set"
- "Search my genome for riboswitches and snoRNAs with Infernal"
- "Do I have enough evidence to call miRNAs, or do I need small-RNA-seq?"

## Example Prompts

### General ncRNA Annotation

> "Run Infernal cmscan against Rfam with GA thresholds and clan competition, then deoverlap, and report counts by family - and remind me this is a lower bound."

> "Annotate ncRNAs in my eukaryotic genome with Infernal plus tRNAscan-SE in eukaryotic mode, and give me the tRNA high-confidence set."

### tRNA-Specific

> "Find tRNA genes and anticodons in my bacterial genome with tRNAscan-SE bacterial mode, separating pseudogenes from the high-confidence set."

> "My organellar tRNAs are being flagged as pseudogenes - run organellar mode."

### Interpretation

> "I have 180 tRNA 'genes' in a bacterium - is that real or fragments/pseudogenes?"

> "Combine Infernal and tRNAscan-SE into one GFF, preferring tRNAscan-SE for tRNAs."

## What the Agent Will Do

1. Press the Rfam CM database (`cmpress`) and record the Rfam release
2. Run Infernal cmscan with `--cut_ga --nohmmonly --fmt 2 --clanin` and `-Z` set from genome size
3. Deoverlap clan hits (`grep -v " = "`) so one locus is not counted several times
4. Run tRNAscan-SE in the correct domain mode and report the high-confidence set
5. Combine results, preferring tRNAscan-SE for tRNA
6. Report counts by class as a floor, flag rDNA copy number as collapsed, and route miRNA/lncRNA to expression evidence

## Tips

- **BLAST is the wrong tool** - Use covariance models (Infernal/Rfam) for structured ncRNA; BLAST is blind to compensatory covariation and misses structure-conserved homologs.
- **`--cut_ga` is correctness, not preference** - Per-family curator thresholds; overriding GA genome-wide to "find more" imports the false positives the curator excluded.
- **Clan competition is mandatory** - `--fmt 2 --clanin` + `grep -v " = "`; without it one locus matching several related models is counted several times.
- **tRNAscan-SE domain mode matters** - `-B`/`-A`/`-E`/`-G`; the wrong mode mis-scores. Report the high-confidence set, not raw hits (SINEs/pseudogenes inflate raw counts 2-10x in eukaryotes).
- **Counts are floors** - rDNA arrays collapse in assemblies (rRNA copy number is off by orders of magnitude); tRNA-gene number is not tRNA abundance.
- **miRNA needs small-RNA-seq** - Genomic hairpin prediction is unreliable; report homology-only hits as candidates. Use plant-specific tools for plants.
- **lncRNA is not a homology task** - It is transcript assembly + coding-potential; never quote a lncRNA count without source + version (RefSeq vs GENCODE vs NONCODE differ by policy, not biology).

## Related Skills

- genome-annotation/prokaryotic-annotation - Bakta/Prokka wrap barrnap + tRNAscan-SE + Infernal
- genome-annotation/eukaryotic-gene-prediction - Protein-coding prediction excludes ncRNAs
- genome-annotation/annotation-qc - tRNA/rRNA count sanity in the QC panel
- rna-structure/ncrna-search - Targeted covariance-model homology searches
