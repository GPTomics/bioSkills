# read-qc

## Overview

Read quality control and preprocessing - the first step in any NGS workflow. Covers raw-FASTQ quality reporting (FastQC/MultiQC), adapter trimming, quality/length filtering, all-in-one fastp preprocessing, contamination screening, UMI processing, and RNA-seq-specific post-alignment QC. Each skill is decision-grade: it teaches the mechanism and the trap behind the choice (why FastQC false-fails on RNA-seq, why a species screen cannot see a same-species swap, why non-UMI RNA-seq must not be deduplicated), not just the tool call.

**Tool type:** cli | **Primary tools:** FastQC, MultiQC, fastp, Cutadapt, Trimmomatic, FastQ Screen, umi_tools, fgbio, RSeQC

## Skills

| Skill | Description |
|-------|-------------|
| quality-reports | Per-file and cross-sample QC with FastQC/falco/MultiQC; read the plots by assay, not the traffic light |
| adapter-trimming | Remove read-through adapter with Cutadapt/Trimmomatic; adapter content is an insert-size readout |
| quality-filtering | Quality/length/N/complexity filtering; trim lightly or not at all, always with a length filter |
| fastp-workflow | All-in-one preprocessing in one pass: overlap adapter trim, filter, poly-G, QC report |
| contamination-screening | Cross-species (FastQ Screen/Kraken2) and same-species swaps (SNP fingerprints) |
| umi-processing | UMI extraction, directional dedup (umi_tools), and single-strand/duplex consensus (fgbio) |
| rnaseq-qc | Post-alignment RNA QC: strandedness, gene-body coverage, TIN, rRNA, read distribution |

## Method Selection

| Decision | Options | Guidance |
|----------|---------|----------|
| Per-file vs cohort QC | FastQC/falco vs MultiQC | FastQC per file; MultiQC makes the cohort the unit of review (outliers, batch effects) |
| Trim adapter vs quality | adapter-trimming vs quality-filtering | Adapter trimming is near-universal; quality trimming is usually unnecessary before soft-clipping aligners and harmful if aggressive |
| One tool vs separate | fastp vs Cutadapt+Trimmomatic+FastQC | fastp for bulk single-pass; Cutadapt for small-RNA/amplicon precision |
| Contamination type | species screen vs SNP fingerprint | Species screen for organisms; SNP fingerprint for same-species swaps/mixtures (orthogonal, need both for human cohorts) |
| Dedup or not | UMI dedup vs coordinate dedup vs none | UMIs for molecule counting; coordinate dedup for non-UMI DNA; NEVER dedup non-UMI bulk RNA-seq |
| RNA-seq strandedness | infer empirically | Always infer (infer_experiment.py / salmon -l A); wrong strand silently halves/zeros counts |

## Example Prompts
- "Run FastQC on all my samples and aggregate with MultiQC"
- "My FastQC fails per-base content in the first 12 bases of RNA-seq, is that a problem?"
- "Remove TruSeq adapters from my paired-end reads and discard reads under 20 bp"
- "Trim small-RNA adapters and keep only 18-30 nt inserts"
- "Run fastp on my NovaSeq data with poly-G trimming"
- "Check my reads for cross-species contamination"
- "Verify my human samples are not swapped or mixed"
- "Extract UMIs and deduplicate after alignment"
- "Build duplex consensus reads for my ctDNA library"
- "Infer the strandedness of my RNA-seq library before quantifying"

## Requirements

```bash
# Quality reporting and aggregation
conda install -c bioconda fastqc multiqc falco seqkit nanoplot

# Trimming and filtering
conda install -c bioconda fastp cutadapt trimmomatic

# Contamination screening
conda install -c bioconda fastq-screen kraken2 bracken bbmap sortmerna somalier

# UMI processing
conda install -c bioconda umi_tools fgbio

# RNA-seq QC
conda install -c bioconda rseqc qualimap rna-seqc picard salmon
```

## Related Skills

- **sequence-io** - FASTQ file reading, writing, and statistics
- **read-alignment** - Downstream mapping (BWA/STAR/Bowtie2/HISAT2)
- **alignment-files** - BAM processing and coordinate-based duplicate handling
- **rna-quantification** - Strand-aware counting after RNA-seq QC
- **metagenomics** - QC before taxonomic classification
- **liquid-biopsy** - Duplex-consensus rare-variant detection
