# STAR RNA-seq Alignment - Usage Guide

## Overview

STAR is the fast, splice-aware RNA-seq aligner that places reads across exon-exon junctions on the genome. It is the choice when the deliverable needs genomic coordinates -- novel isoforms, fusions, RNA variant calling, coverage tracks, splicing QC, or single-cell counting (STARsolo) -- and has memory to spare (~30 GB for human). The decisions that actually shape the result are the splice-junction database (built from a GTF at sjdbOverhang = readlength - 1), the per-sample vs cohort two-pass choice, the 255-for-unique MAPQ that breaks GATK, and the strand column STAR reports; this skill emphasizes those over the basic command.

## Prerequisites

```bash
conda install -c bioconda star samtools
```

- A genome FASTA and a matching GTF annotation from the same provider/release (contig naming must agree; reconcile with alignment-files).
- ~30 GB RAM for a human index; route memory-constrained jobs to hisat2-alignment.
- The read length, to set `--sjdbOverhang` at index build.

## Quick Start

Tell your AI agent what you want to do:
- "Build a STAR index for my genome and GTF with the right sjdbOverhang for 150 bp reads"
- "Align my paired-end RNA-seq reads with STAR in two-pass mode and get gene counts"
- "Detect my library strandedness from the STAR GeneCounts output"
- "Align RNA-seq for GATK variant calling without the empty-VCF MAPQ problem"
- "Run STAR for fusion detection with STAR-Fusion"

## Example Prompts

### Index generation
> "Create a STAR index for genome.fa with genes.gtf for 150 bp reads, and reduce genomeSAindexNbases appropriately because my genome is a 5 Mb bacterial assembly."

### Two-pass and counts
> "Align my paired-end RNA-seq with STAR two-pass, output a coordinate-sorted BAM and gene counts, and tell me from the counts whether my library is stranded."

### RNA variant calling
> "Align RNA-seq with STAR for GATK variant calling and make sure the unique-read MAPQ does not cause an empty VCF."

### Cohort splicing
> "I am doing differential splicing across 30 samples -- how do I run STAR two-pass so the junction set is identical across samples instead of per-sample?"

### Fusion detection
> "Run STAR with chimeric output for STAR-Fusion on my tumor RNA-seq."

## What the Agent Will Do

1. Build the index with `--sjdbGTFfile` and `--sjdbOverhang = readlength - 1`, reducing `--genomeSAindexNbases` for small genomes.
2. Align to a coordinate-sorted BAM (STAR sorts natively; a later samtools sort is redundant), adding `--quantMode GeneCounts` to also read strandedness.
3. Choose per-sample two-pass for novel-junction work, or the cohort pooling recipe (shared SJ.out.tab) for splicing comparisons.
4. Set `--outSAMmapqUnique 60` for any STAR -> GATK path so the 255 uniques are not dropped.
5. Avoid MAPQ-filtering RNA output for counting (it deletes multimappers), and route multimapper handling to the counter.
6. Reconcile genome/GTF contig naming and run the QC gate (alignment-files/bam-statistics) before counting or calling.

## Index sjdbOverhang by Read Length

| Read length | sjdbOverhang |
|-------------|--------------|
| 50 bp | 49 |
| 75 bp | 74 |
| 100 bp | 99 |
| 150 bp | 149 |

## Strandedness from GeneCounts

| Library type | ReadsPerGene column | Equivalent counter setting |
|--------------|---------------------|----------------------------|
| Unstranded | column 2 | htseq `-s no` / featureCounts `-s 0` |
| Forward (Ligation) | column 3 | htseq `-s yes` / featureCounts `-s 1` |
| Reverse (dUTP, TruSeq -- the common case) | column 4 | htseq `-s reverse` / featureCounts `-s 2` |

## Tips

- Set `--sjdbOverhang` to read length - 1 at index build; the default 100 quietly degrades junction sensitivity for short reads.
- STAR already coordinate-sorts with `--outSAMtype BAM SortedByCoordinate` -- do not run `samtools sort` again.
- For any STAR -> GATK RNA-variant path, set `--outSAMmapqUnique 60` so the 255 uniques are not dropped as "unavailable," and inject read groups with `--outSAMattrRGline ID:.. SM:.. PL:ILLUMINA LB:..` (space-separated tags, not bwa's `-R '@RG\t...'`) because GATK requires them.
- STAR does not auto-detect gzip; pass `--readFilesCommand zcat` for `.gz` input (bwa/bowtie2/HISAT2 auto-detect it). For htseq-count, emit `--outSAMtype BAM Unsorted` or re-sort with `samtools sort -n`; featureCounts accepts coordinate order.
- Do not MAPQ-filter STAR output for counting; a `-q 10` filter deletes every multimapper and under-counts gene families.
- For cohort splicing/sQTL, pool the pass-1 SJ.out.tab files and re-align all samples against one common junction set; per-sample two-pass is a batch effect.
- Infer strandedness from the GeneCounts columns rather than assuming; the wrong column roughly halves the counts.
- Reduce `--genomeSAindexNbases` for small genomes or STAR silently builds a broken index.

## STAR vs HISAT2

| Factor | STAR | HISAT2 |
|--------|------|--------|
| Memory (human) | ~30 GB | ~7 GB |
| Speed | very fast | fast |
| Native gene counts | yes (GeneCounts) | no |
| Fusion detection | yes (chimeric) | no |
| MAPQ for GATK | 255 -> set 60 | 60 (GATK-friendly) |
| Downstream coupling | featureCounts, RSEM, STAR-Fusion, STARsolo | StringTie/Cufflinks via --dta |

## Related Skills

- hisat2-alignment - Low-memory splice-aware alternative to STAR
- bwa-alignment - DNA short-read mapping (when reads do not cross junctions)
- read-qc/rnaseq-qc - RNA destination metrics: rRNA, gene-body coverage, strandedness
- read-qc/fastp-workflow - Trim adapters/poly-A before alignment
- alignment-files/bam-statistics - flagstat/idxstats QC gate; what a high mapping rate hides; contig naming
- rna-quantification/featurecounts-counting - Count aligned reads over genes (NH-aware multimapper handling)
- rna-quantification/alignment-free-quant - Salmon/kallisto when only known-transcript DE is needed
- differential-expression/deseq2-basics - Downstream DE from the count matrix
- single-cell/data-io - STARsolo single-cell counts into a single-cell workflow
