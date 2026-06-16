# HISAT2 RNA-seq Alignment - Usage Guide

## Overview

HISAT2 is a memory-efficient splice-aware RNA-seq aligner. Its hierarchical graph FM-index does spliced alignment at roughly a quarter of STAR's memory (~7 GB vs ~30 GB for human), its MAPQ is GATK-friendly (60 for unique, no 255 problem), and its SNP/haplotype graph index can reduce reference bias before any read is mapped. It is the right choice on a memory-constrained machine and the standard front end for StringTie/Cufflinks transcript assembly (via --dta). This skill emphasizes the memory and graph-index advantages, the --dta trade-off, and the strandedness setting that silently halves counts when wrong.

## Prerequisites

```bash
conda install -c bioconda hisat2 samtools
```

- A HISAT2 index (built with `hisat2-build`) or a prebuilt index (`grch38_tran`, `grch38_snp_tran`) from the HISAT2 site; a full human annotation-aware build needs substantial RAM.
- The library strandedness (or a plan to infer it); set `--rna-strandness` accordingly.

## Quick Start

Tell your AI agent what you want to do:
- "Align my RNA-seq reads with HISAT2 on a low-memory machine"
- "Build a HISAT2 index with splice sites and exons from my GTF"
- "Run HISAT2 with the right strandedness for my TruSeq library"
- "Align for StringTie transcript assembly with the --dta flag"
- "Use the SNP-aware graph index to reduce reference bias"

## Example Prompts

### Index building
> "Build a HISAT2 index for genome.fa with splice sites and exons from genes.gtf, and tell me if I should use a prebuilt index instead to avoid the RAM cost."

### Basic alignment
> "Align my paired-end TruSeq stranded RNA-seq with HISAT2 using the correct --rna-strandness, and output a sorted, indexed BAM."

### Transcript assembly
> "Align reads for StringTie with HISAT2 --dta, and explain what --dta trades away for plain counting."

### Memory choice
> "STAR needs ~30 GB and my node has 16 GB -- align my RNA-seq with HISAT2 instead and explain the trade-off."

## What the Agent Will Do

1. Build or locate the index (plain, annotation-aware, or SNP-graph), preferring a prebuilt index when a full human `--ss --exon` build would exhaust RAM.
2. Set `--rna-strandness` from the library (RF for the common dUTP/TruSeq case), inferring it if unknown.
3. Add `--dta` only when StringTie/Cufflinks assembly is the downstream, never for plain counting.
4. Add read groups and stream to a coordinate-sorted, indexed BAM.
5. For novel-junction work, run a manual two-pass and merge junctions across a cohort to avoid a per-sample batch effect.
6. Reconcile genome/GTF contig naming and run the QC gate (alignment-files/bam-statistics) before counting.

## Strandedness Settings

| Kit / method | HISAT2 flag |
|--------------|-------------|
| Unstranded | (omit --rna-strandness) |
| Illumina TruSeq Stranded mRNA | --rna-strandness RF |
| dUTP method | --rna-strandness RF |
| Ligation / forward method | --rna-strandness FR |
| Single-end reverse / forward | --rna-strandness R / F |

## Tips

- Reach for HISAT2 when memory is the constraint; its ~7 GB graph index runs where STAR's ~30 GB does not.
- Set `--rna-strandness` (RF for most Illumina stranded libraries) or sense reads fall in "no feature" and counts roughly halve; infer it with RSeQC infer_experiment.py if unknown.
- Use `--dta` only for transcript assembly (StringTie/Cufflinks); it discards short-anchor junction reads, hurting plain counting.
- Use a prebuilt `grch38_tran` / `grch38_snp_tran` index, or pass junctions at align time via `--known-splicesite-infile`, to avoid the heavy RAM of a full annotation-aware build.
- The SNP-graph index removes reference bias for known variants in the index itself; private/novel-variant ASE still needs WASP or a personalized reference.
- HISAT2 gives unique reads MAPQ 60, so its output goes straight into GATK without the MAPQ reassignment STAR requires.

## HISAT2 vs STAR

| Factor | Choose HISAT2 | Choose STAR |
|--------|---------------|-------------|
| Memory | limited (<32 GB) | available (>=32 GB) |
| Native gene counts | not needed | needed (GeneCounts) |
| Fusion detection | not needed | needed (chimeric) |
| Transcript assembly | StringTie/Cufflinks via --dta | -- |
| MAPQ for GATK | 60 (friendly) | 255 -> set 60 |
| Novel-junction sensitivity | good | highest (2-pass) |

## Related Skills

- star-alignment - Feature-rich, higher-RAM splice-aware alternative (native counts, fusions)
- bwa-alignment - DNA short-read mapping (when reads do not cross junctions)
- read-qc/rnaseq-qc - RNA destination metrics: rRNA, gene-body coverage, strandedness
- read-qc/fastp-workflow - Trim adapters/poly-A before alignment
- alignment-files/bam-statistics - flagstat/idxstats QC gate; what a high mapping rate hides; contig naming
- rna-quantification/featurecounts-counting - Count aligned reads over genes
- rna-quantification/alignment-free-quant - Salmon/kallisto when only known-transcript DE is needed
- differential-expression/deseq2-basics - Downstream DE from the count matrix
