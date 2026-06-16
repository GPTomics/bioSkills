# read-alignment

## Overview

Map short reads to a reference genome with the right aligner for the assay, on a correctly-chosen reference: DNA alignment with bwa-mem2/Bowtie2 and RNA-seq spliced alignment with STAR/HISAT2. Each skill carries its own reference, MAPQ, and read-group decisions; the cross-tool MAPQ scale, contig-naming reconciliation, and the BAM QC gate are owned by alignment-files.

**Tool type:** cli | **Primary tools:** bwa-mem2, bowtie2, STAR, HISAT2, samtools

## Skills

| Skill | Description |
|-------|-------------|
| bwa-alignment | DNA short reads with bwa-mem2: read groups, ALT/decoy-aware mapping, the dedup ordering, -M/-Y, -K determinism |
| bowtie2-alignment | ChIP/ATAC/CUT&RUN with Bowtie2: end-to-end vs local, sensitivity presets, fragment-geometry flags |
| star-alignment | RNA spliced alignment with STAR: sjdbOverhang, two-pass, GeneCounts strandedness, the 255 MAPQ, fusions |
| hisat2-alignment | Low-memory RNA spliced alignment with HISAT2: graph FM-index, SNP-graph, strandedness, --dta |

## Method Selection

| If you need to... | Use | Why |
|-------------------|-----|-----|
| Align DNA (WGS/WES/somatic) for variant calling | bwa-alignment | seed-and-extend standard; ALT/decoy-aware; GATK/DeepVariant expect its output |
| Align ChIP-seq / ATAC-seq / CUT&RUN | bowtie2-alignment | end-to-end vs local plus fragment-geometry flags feed the peak caller |
| Align RNA-seq with ample RAM / want fusions or native counts | star-alignment | fastest splice-aware; GeneCounts, chimeric output, two-pass |
| Align RNA-seq on limited RAM or for StringTie | hisat2-alignment | ~7 GB graph FM-index; --dta for transcript assembly; GATK-friendly MAPQ |
| Get a DE count matrix on known transcripts without a BAM | route to rna-quantification/alignment-free-quant | Salmon/kallisto are faster and model multimapping; do not align |
| Interpret BAM QC (mapping rate, dup, contig naming), or diagnose zero counts / empty VCF | route to alignment-files | flagstat/idxstats interpretation, the cross-tool MAPQ scale, and reference-naming reconciliation live there |

## Example Prompts

- "Which aligner should I use for my ChIP-seq / RNA-seq / WGS data and why?"
- "Pick the right GRCh38 analysis set (decoy/ALT) for human variant calling"
- "My featureCounts output is all zeros despite a 95% mapping rate - what is wrong?"
- "Align my paired-end WGS reads with bwa-mem2 and read groups for GATK"
- "Map reads for structural-variant calling with soft-clipped supplementary alignments"
- "Align ChIP-seq reads with Bowtie2 and filter multimappers"
- "Align ATAC-seq reads with the right fragment and soft-clipping settings"
- "Align RNA-seq with STAR in two-pass mode and detect my library strandedness"
- "Align RNA-seq for GATK variant calling without the empty-VCF MAPQ problem"
- "Align RNA-seq with HISAT2 on a 16 GB machine and explain the STAR trade-off"

## Requirements

```bash
conda install -c bioconda bwa-mem2 bwa bowtie2 star hisat2 samtools
```

## Related Skills

- **read-qc** - Upstream quality control and trimming before alignment
- **alignment-files** - Post-alignment BAM sort/dedup/index/stats mechanics
- **rna-quantification** - Downstream counting / alignment-free quantification of aligned reads
- **variant-calling** - Variant calling from aligned reads
- **chip-seq** / **atac-seq** - Peak calling from ChIP/ATAC alignments
- **methylation-analysis** - Bisulfite alignment (Bismark wraps Bowtie2)
